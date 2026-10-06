#!/usr/bin/env python3
"""Station coordinate adjustment from the fixed network residuals of nrtk-net.

A station position error dr enters every ambiguity-fixed double-difference
residual of its baselines as -(e_s - e_r) . dr (e: unit vectors station ->
satellite). With one station held fixed, a least-squares fit of the
non-dispersive residuals of a whole day gives the other stations' corrections
to a few mm, using nrtk's own integer fixing instead of a rover-style filter.
The remaining relative wet delay of each baseline (the residuals carry the
Saastamoinen standard-atmosphere model only, not the filter's wet-delay state)
is co-estimated as one value per hour and eliminated group by group (Schur
complement), so the normal equations only hold the coordinates: a 120-station
network with 1000 baselines takes a few minutes. Wrong fixes (residuals off by
whole wavelengths) are removed by iterative reweighting. The elevation weights
sin^2(el) are relative to the robust residual sigma, and the priors on the wet
delay (--zwd-prior) and the coordinates are scaled with that sigma; a prior in
absolute units against unit-variance weights would shrink the hourly wet delays
towards zero and let the unabsorbed part leak into the heights (2-3 cm on the
Taiwan network before this was fixed).

  nrtk-net -c conf.yaml -o work --sat-csv            # residual CSVs blsat_A-B.csv
  tools/twdata/adjust_coords.py -c conf.yaml -r work --fix LNK2 [--apply conf_new.yaml]

Repeat (nrtk-net with the corrected coordinates) until the corrections are
within a millimetre; a second day processed with the same coordinates should
then also give corrections at the millimetre level (vertical: ~1 cm).
"""
import argparse, csv, glob, json, math, os, re, subprocess, sys
from multiprocessing import Pool
import numpy as np

GRAD_C=0.0032   # Chen & Herring gradient mapping constant

def confjson(bindir, conf):
    return json.loads(subprocess.run([os.path.join(bindir,'nrtk-tools'),'confjson',conf],
                                     capture_output=True,text=True,check=True).stdout)

def unit(az,el):
    return (math.cos(el)*math.sin(az), math.cos(el)*math.cos(az), math.sin(el))

def mgrad(el):
    return 1.0/(math.sin(el)*math.tan(el)+GRAD_C)

def parse_file(args):
    """one blsat CSV -> arrays: d (n,3) unit-vector difference, m (n,3) troposphere partials
    (relative ZWD, north and east gradient), y (n), w (n), hour (n)"""
    fn,hours,elmin,w0,w1=args
    ep={}
    for r in csv.DictReader(open(fn)):
        if r['fixed']!='1' or r['newfix']=='1': continue
        ep.setdefault(r['time'],[]).append(r)
    D=[];M=[];Y=[];W=[];H=[];S=[]
    for t,sats in ep.items():
        refs={r['sat'][0]:r for r in sats if r['isref']=='1'}
        if not refs: continue
        hh,mm,ss=[float(x) for x in t[11:].split(':')]
        th=hh+mm/60+ss/3600
        if th<w0 or th>=w1: continue
        hour=int(th//hours)
        for r in sats:
            if r['isref']=='1' or r['sat'][0] not in refs: continue
            rr=refs[r['sat'][0]]
            el=float(r['el_m']); elr=float(rr['el_m'])
            if el<elmin: continue
            azs=math.radians(float(r['az_m'])); azr=math.radians(float(rr['az_m']))
            els=math.radians(el); elr_=math.radians(elr)
            es=unit(azs,els); er=unit(azr,elr_)
            D.append((es[0]-er[0],es[1]-er[1],es[2]-er[2]))
            gs=mgrad(els); gr=mgrad(elr_)
            M.append((1.0/math.sin(els)-1.0/math.sin(elr_),
                      gs*math.cos(azs)-gr*math.cos(azr), gs*math.sin(azs)-gr*math.sin(azr)))
            Y.append(float(r['dd_nondisp'])-float(r['dd_trop_model']))
            W.append(math.sin(math.radians(el))**2)
            H.append(hour); S.append(r['sat'])
    return (np.array(D,float).reshape(-1,3),np.array(M,float).reshape(-1,3),np.array(Y,float),np.array(W,float),np.array(H,int),np.array(S))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('-c','--conf',required=True)
    ap.add_argument('-r','--results',required=True,help='nrtk-net output dir with blsat_*.csv')
    ap.add_argument('-b','--bin',default='build/bin')
    ap.add_argument('--fix',required=True,help='station held fixed (datum)')
    ap.add_argument('--hours',type=float,default=1.0,help='relative wet delay piece length (h)')
    ap.add_argument('--elmin',type=float,default=15.0,help='elevation cutoff (deg)')
    ap.add_argument('--gradients',action='store_true',help='co-estimate relative tropospheric gradients per piece')
    ap.add_argument('--zwd-prior',type=float,default=0.05,help='prior std of the relative ZWD per piece (m)')
    ap.add_argument('--grad-prior',type=float,default=0.005,help='prior std of the relative gradients per piece (m)')
    ap.add_argument('--window',type=float,nargs=2,default=[0.0,24.0],metavar=('H0','H1'),help='use epochs in [H0,H1) hours of day')
    ap.add_argument('--satstats',action='store_true',help='print per-satellite statistics of the final residuals')
    ap.add_argument('--apply',help='write the configuration with corrected coordinates here')
    ap.add_argument('-j','--jobs',type=int,default=8)
    a=ap.parse_args()
    cfg=confjson(a.bin,a.conf)
    sta={s['id']:s for s in cfg['stations']}
    if a.fix not in sta: sys.exit('unknown station '+a.fix)
    files=sorted(glob.glob(os.path.join(a.results,'blsat_*.csv')))
    bls=[]
    for fn in files:
        name=os.path.basename(fn)[6:-4]
        m,aa=name.split('-')
        if m in sta and aa in sta: bls.append((fn,m,aa))
    if not bls: sys.exit('no blsat_*.csv in '+a.results)
    free=sorted({s for _,m,aa in bls for s in (m,aa)}-{a.fix})
    col={s:3*i for i,s in enumerate(free)}
    nsta=3*len(free)
    with Pool(a.jobs) as pool:
        parsed=pool.map(parse_file,[(fn,a.hours,a.elmin,a.window[0],a.window[1]) for fn,_,_ in bls])
    ntot=sum(len(p[2]) for p in parsed)
    print("%d residuals from %d baselines, %d stations free"%(ntot,len(bls),len(free)),file=sys.stderr)

    # groups: (baseline, hour) with coordinate design A (n x 6: aux -d, master +d)
    nt=3 if a.gradients else 1
    groups=[]
    for (fn,m,aa),(D,M,Y,W,H,S) in zip(bls,parsed):
        if len(Y)==0: continue
        ca=col.get(aa,-1); cm=col.get(m,-1)
        for h in np.unique(H):
            sel=H==h
            groups.append((ca,cm,D[sel],M[sel][:,:nt],Y[sel],W[sel],S[sel]))
    x=np.zeros(nsta); s0=0.02
    wscale=[np.ones(len(g[4])) for g in groups]
    for it in range(4):
        # priors relative to the residual sigma (the weights W are unitless sin^2(el))
        prior_t=np.diag([(s0/a.zwd_prior)**2]+[(s0/a.grad_prior)**2]*(nt-1))
        N=np.eye(nsta)*(s0/1.0)**2   # 1 m prior keeps stations without data finite
        b=np.zeros(nsta)
        for g,ws in zip(groups,wscale):
            ca,cm,D,M,Y,W,_=g
            w=W*ws
            cols=[];A=[]
            if ca>=0: cols+=[ca,ca+1,ca+2]; A.append(-D)
            if cm>=0: cols+=[cm,cm+1,cm+2]; A.append(D)
            if not cols: continue
            A=np.hstack(A)                        # n x k
            AtW=A.T*w                             # k x n
            Ncc=AtW@A; bc=AtW@Y
            Ntt=(M.T*w)@M+prior_t                 # nt x nt
            Nct=AtW@M; bt=(M.T*w)@Y               # k x nt, nt
            Nti=np.linalg.inv(Ntt)
            Ncc-=Nct@Nti@Nct.T; bc-=Nct@(Nti@bt)
            N[np.ix_(cols,cols)]+=Ncc; b[cols]+=bc
        x=np.linalg.solve(N,b)
        # residuals with the eliminated troposphere, robust reweighting
        allres=[]
        for g,ws in zip(groups,wscale):
            ca,cm,D,M,Y,W,_=g
            w=W*ws
            Ax=np.zeros(len(Y))
            if ca>=0: Ax+=-D@x[ca:ca+3]
            if cm>=0: Ax+=D@x[cm:cm+3]
            t=np.linalg.solve((M.T*w)@M+prior_t,(M.T*w)@(Y-Ax))
            allres.append(Y-Ax-M@t)
        res=np.concatenate(allres)
        s0=1.4826*np.median(np.abs(res))
        rej=0
        for gi,r in enumerate(allres):
            ok=np.abs(r)<3*s0
            wscale[gi]=ok.astype(float); rej+=int((~ok).sum())
        print("iteration %d: robust sigma %.1f mm, %d residuals rejected"%(it,1e3*s0,rej),file=sys.stderr)
    if a.satstats:
        st={}
        for g,r,ws in zip(groups,allres,wscale):
            for sat in np.unique(g[6]):
                sel=g[6]==sat
                d=st.setdefault(sat,[0,0.0,0.0,0])
                rr=r[sel]; ok=ws[sel]>0
                d[0]+=int(ok.sum()); d[1]+=float(rr[ok].sum()); d[2]+=float((rr[ok]**2).sum()); d[3]+=int((~ok).sum())
        print("satellite statistics of the final residuals (mm): sat n mean rms rejected%",file=sys.stderr)
        for sat,(n,s1,s2,nrej) in sorted(st.items(),key=lambda kv:-abs(kv[1][1]/max(1,kv[1][0]))):
            if n: print("  %s %8d %+6.1f %6.1f %5.1f%%"%(sat,n,1e3*s1/n,1e3*math.sqrt(s2/n),100*nrej/(n+nrej)),file=sys.stderr)
    cov=np.linalg.inv(N)*s0**2
    print("station  dE    dN    dU   (mm, +-formal)  [%s held fixed]"%a.fix)
    corr={}
    for s in free:
        c=col[s]; d=x[c:c+3]*1e3; sd=np.sqrt(np.diag(cov)[c:c+3])*1e3
        print("%s %6.1f %6.1f %6.1f   +- %4.1f %4.1f %4.1f"%(s,d[0],d[1],d[2],sd[0],sd[1],sd[2]))
        corr[s]=x[c:c+3]
    if a.apply:
        txt=open(a.conf).read()
        for s,d in corr.items():
            x0,y0,z0=sta[s]['ecef']; lat,lon=[math.radians(v) for v in sta[s]['llh'][:2]]
            sl,cl,sp,cp=math.sin(lon),math.cos(lon),math.sin(lat),math.cos(lat)
            e,n,u=d
            dx=-sl*e-sp*cl*n+cp*cl*u; dy=cl*e-sp*sl*n+cp*sl*u; dz=cp*n+sp*u
            mm=re.search(r'(- id: %s\b[^\n]*\n(?:\s*enabled:[^\n]*\n)?\s*ecef: \[)([^\]]*)\]'%re.escape(s),txt)
            if not mm: print("# %s: not found in config"%s,file=sys.stderr); continue
            txt=txt[:mm.start(2)]+"%.4f, %.4f, %.4f"%(x0+dx,y0+dy,z0+dz)+txt[mm.end(2):]
        open(a.apply,'w').write(txt)
        print("wrote",a.apply,file=sys.stderr)
if __name__=='__main__':
    main()
