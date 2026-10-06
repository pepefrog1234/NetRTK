#!/usr/bin/env python3
"""Station coordinate adjustment from the fixed network residuals of nrtk-net.

A station position error dr enters every ambiguity-fixed double-difference
residual of its baselines as -(e_s - e_r) . dr (e: unit vectors station ->
satellite). With one station held fixed, a least-squares fit of the
non-dispersive residuals of a whole day (tens of thousands per baseline) gives
the other stations' corrections to a few mm, using nrtk's own integer fixing
instead of a rover-style filter. The remaining relative wet delay of each
baseline is co-estimated as one value per hour. Wrong fixes (residuals off by
whole wavelengths) are removed by iterative reweighting.

  nrtk-net -c conf.yaml -o work --sat-csv            # residual CSVs blsat_A-B.csv
  tools/twdata/adjust_coords.py -c conf.yaml -r work --fix LNK2 [--apply conf_new.yaml]

Repeat (nrtk-net with the corrected coordinates) until the corrections are
within their formal errors; a second day processed with the same coordinates
should then also give corrections near zero.
"""
import argparse, csv, glob, json, math, os, subprocess, sys
import numpy as np

def confjson(bindir, conf):
    return json.loads(subprocess.run([os.path.join(bindir,'nrtk-tools'),'confjson',conf],
                                     capture_output=True,text=True,check=True).stdout)

def unit(az,el):
    return np.array([math.cos(el)*math.sin(az), math.cos(el)*math.cos(az), math.sin(el)])

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('-c','--conf',required=True)
    ap.add_argument('-r','--results',required=True,help='nrtk-net output dir with blsat_*.csv')
    ap.add_argument('-b','--bin',default='build/bin')
    ap.add_argument('--fix',required=True,help='station held fixed (datum)')
    ap.add_argument('--elmin',type=float,default=15.0)
    ap.add_argument('--hours',type=float,default=1.0,help='relative wet delay piece length (h)')
    ap.add_argument('--apply',help='write the configuration with corrected coordinates here')
    a=ap.parse_args()
    cfg=confjson(a.bin,a.conf)
    sta={s['id']:s for s in cfg['stations']}
    ids=[s['id'] for s in cfg['stations'] if s['enabled']]
    if a.fix not in sta: sys.exit('unknown station '+a.fix)
    free=[s for s in ids if s!=a.fix]
    col={s:3*i for i,s in enumerate(free)}
    nsta=3*len(free)

    rows=[]   # (coef dict, value, weight)
    files=sorted(glob.glob(os.path.join(a.results,'blsat_*.csv')))
    if not files: sys.exit('no blsat_*.csv in '+a.results)
    ntrop=0; blidx={}
    for fn in files:
        name=os.path.basename(fn)[6:-4]; m,aa=name.split('-')
        if m not in sta or aa not in sta: continue
        blidx[name]=ntrop; ntrop+=1
    def tcol(name,hour): return nsta+blidx[name]*int(round(24/a.hours))+int(hour//a.hours)
    ncol=nsta+ntrop*int(round(24/a.hours))
    N=np.zeros((ncol,ncol)); b=np.zeros(ncol)
    obs=[]
    for fn in files:
        name=os.path.basename(fn)[6:-4]; m,aa=name.split('-')
        if name not in blidx: continue
        ep={}
        for r in csv.DictReader(open(fn)):
            if r['fixed']!='1' or r['newfix']=='1': continue
            ep.setdefault(r['time'],[]).append(r)
        for t,sats in ep.items():
            # every satellite refers to the reference of its own system (letter of the id)
            refs={r['sat'][0]:r for r in sats if r['isref']=='1'}
            if not refs: continue
            hh,mm,ss=[float(x) for x in t[11:].split(':')]
            hour=hh+mm/60+ss/3600
            for r in sats:
                if r['isref']=='1' or r['sat'][0] not in refs: continue
                rr=refs[r['sat'][0]]
                el_a=float(r['el_m'])  # elevation at master side (close to the aux side at 50 km)
                if el_a<a.elmin: continue
                es=unit(math.radians(float(r['az_m'])),math.radians(el_a))
                er=unit(math.radians(float(rr['az_m'])),math.radians(float(rr['el_m'])))
                d=es-er
                # residual = dd_nondisp - dd_trop_model = -(d.dr_a) + (d.dr_m) + dzwd*(mw_s-mw_r)
                y=float(r['dd_nondisp'])-float(r['dd_trop_model'])
                mw=1.0/math.sin(math.radians(el_a))-1.0/math.sin(math.radians(float(rr['el_m'])))
                coef={}
                if aa in col:
                    for k in range(3): coef[col[aa]+k]=-d[k]
                if m in col:
                    for k in range(3): coef[col[m]+k]=d[k]
                coef[tcol(name,hour)]=mw
                w=math.sin(math.radians(el_a))**2
                obs.append((coef,y,w))
    print("%d residuals from %d baselines, %d stations free"%(len(obs),ntrop,len(free)),file=sys.stderr)
    x=np.zeros(ncol); wts=np.array([o[2] for o in obs])
    for it in range(4):
        N[:]=0; b[:]=0
        for (coef,y,w),wi in zip(obs,wts):
            idx=list(coef); c=np.array([coef[i] for i in idx])
            N[np.ix_(idx,idx)]+=wi*np.outer(c,c); b[idx]+=wi*c*y
        # weak priors keep the system regular: 5 cm on the tropospheric pieces
        # (empty hours), 1 m on the coordinates (a station without any fixed
        # residual this day is reported with +-1000 mm)
        for i in range(nsta,ncol): N[i,i]+=1.0/0.05**2
        for i in range(nsta): N[i,i]+=1.0
        x=np.linalg.solve(N,b)
        res=np.array([y-sum(coef[i]*x[i] for i in coef) for coef,y,w in obs])
        s0=1.4826*np.median(np.abs(res))
        wts=np.array([o[2] for o in obs])*np.where(np.abs(res)<3*s0,1.0,0.0)
        print("iteration %d: robust sigma %.1f mm, %d residuals rejected"%(it,1e3*s0,int((wts==0).sum())),file=sys.stderr)
    cov=np.linalg.inv(N)*s0**2
    print("station  dE    dN    dU   (mm, +-formal)  [%s held fixed]"%a.fix)
    corr={}
    for s in free:
        c=col[s]; d=x[c:c+3]*1e3; sd=np.sqrt(np.diag(cov)[c:c+3])*1e3
        print("%s %6.1f %6.1f %6.1f   +- %4.1f %4.1f %4.1f"%(s,d[0],d[1],d[2],sd[0],sd[1],sd[2]))
        corr[s]=x[c:c+3]
    if a.apply:
        import re
        txt=open(a.conf).read()
        for s,d in corr.items():
            x0,y0,z0=sta[s]['ecef']; lat,lon,_=[math.radians(v) for v in sta[s]['llh'][:2]]+[0]
            sl,cl,sp,cp=math.sin(lon),math.cos(lon),math.sin(lat),math.cos(lat)
            e,n,u=d
            dx=-sl*e-sp*cl*n+cp*cl*u; dy=cl*e-sp*sl*n+cp*sl*u; dz=cp*n+sp*u
            pat=re.compile(r'(- id: %s\n\s*ecef: \[)([^\]]*)\]'%s)
            mm=pat.search(txt)
            if not mm: print("# %s: not found in config"%s,file=sys.stderr); continue
            txt=txt[:mm.start(2)]+"%.4f, %.4f, %.4f"%(x0+dx,y0+dy,z0+dz)+txt[mm.end(2):]
        open(a.apply,'w').write(txt)
        print("wrote",a.apply,file=sys.stderr)
if __name__=='__main__':
    main()
