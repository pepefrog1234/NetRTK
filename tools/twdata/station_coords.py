#!/usr/bin/env python3
"""Relative station coordinates for a network from daily RINEX files (RTKLIB rnx2rtkp).

The CWA GDMS archive carries only metre-level APPROX POSITION headers. Network
RTK needs the stations' *relative* positions to a few mm; the absolute datum
matters far less (a common shift moves the VRS and the leave-out rover's truth
together). Every station is solved as a baseline from one seed station over a
full day (broadcast ephemerides, troposphere estimated, antenna PCO/PCV from
ANTEX with the types of the RINEX headers, solid tides on):
  static-iflc (default): 24 h static ionosphere-free L1+L2 float solution,
      forward+backward combined; no ambiguity fixing, so the daytime ionosphere
      of a low-latitude 20-50 km baseline cannot produce wrong fixes.
  kinematic-fix: median of the fixed kinematic epochs (30 s data at low
      latitude: 1-7 cm day-to-day scatter, kept for comparison).
Run it for several days and merge with --merge to check day-to-day
consistency. The seed keeps its header position: its error is a common offset
of the whole network.

  station_coords.py -b build/bin --atx igs20.atx --nav BRDC_DDD.rnx --seed SEED.26o.gz
                    -o work/DDD STATION_FILES...          -> work/DDD/coords.csv
  station_coords.py --merge work/270/coords.csv work/271/coords.csv ... [--yaml]

Baselines longer than ~25 km fix poorly on 30 s data in the daytime ionosphere;
solve such stations in a second step from an already solved nearer station
(--seed that station's file, --seed-xyz its solved position).
"""
import argparse, csv, gzip, math, os, statistics, subprocess, sys

def header(fn):
    op=gzip.open if fn.endswith('.gz') else open
    h={'name':None,'xyz':None,'ant':'','rec':'','hgt':0.0,'date':None}
    with op(fn,'rt',errors='replace') as f:
        for line in f:
            lab=line[60:].strip()
            if lab=='MARKER NAME': h['name']=line[:60].strip()
            elif lab=='APPROX POSITION XYZ': h['xyz']=[float(line[i:i+14]) for i in (0,14,28)]
            elif lab=='ANT # / TYPE': h['ant']=' '.join(line[20:40].split())
            elif lab=='REC # / TYPE / VERS': h['rec']=line[20:40].strip()
            elif lab=='ANTENNA: DELTA H/E/N': h['hgt']=float(line[:14])
            elif lab=='TIME OF FIRST OBS': h['date']=[int(float(x)) for x in line[:18].split()]
            elif lab=='END OF HEADER': break
    return h

def write_conf(path, atx, seed_xyz, elmask, navsys, mode):
    static=mode.startswith('static')
    opts={'pos1-posmode':'static' if static else 'kinematic','pos1-frequency':'l1+l2',
          'pos1-soltype':'combined' if static else 'forward',
          'pos1-elmask':str(elmask),'pos1-dynamics':'off','pos1-tidecorr':'on',
          'pos1-ionoopt':'iflc' if mode=='static-iflc' else 'off','pos1-tropopt':'est-ztd','pos1-sateph':'brdc',
          'pos1-navsys':str(navsys),
          'pos2-armode':'off' if mode=='static-iflc' else ('fix-and-hold' if mode=='static-fix' else 'continuous'),
          'pos2-arthres':'3','pos2-elmaskar':'15',
          'pos2-arfilter':'on',
          'pos2-maxage':'30','out-solformat':'xyz','out-outhead':'on','out-timeform':'tow',
          'out-timesys':'gpst','out-outstat':'off',
          'ant1-postype':'single','ant1-anttype':'*','ant2-postype':'xyz',
          'ant2-pos1':'%.4f'%seed_xyz[0],'ant2-pos2':'%.4f'%seed_xyz[1],'ant2-pos3':'%.4f'%seed_xyz[2],
          'ant2-anttype':'*','file-rcvantfile':atx}
    with open(path,'w') as f:
        for k,v in opts.items(): f.write('%-20s=%s\n'%(k,v))

def read_pos(fn):
    out=[]
    for l in open(fn):
        if l.startswith('%') or not l.strip(): continue
        p=l.split()
        if len(p)<7: continue
        out.append((float(p[1]),float(p[2]),float(p[3]),float(p[4]),int(p[5]),int(p[6])))
    return out

def solve(a):
    os.makedirs(a.out,exist_ok=True)
    seed=header(a.seed)
    if a.seed_xyz: seed['xyz']=list(a.seed_xyz)
    conf=os.path.join(a.out,'rtk.conf')
    write_conf(conf,a.atx,seed['xyz'],a.elmask,a.navsys,a.mode)
    rows=[]
    for fn in a.files:
        h=header(fn)
        if os.path.abspath(fn)==os.path.abspath(a.seed):
            rows.append(dict(id=h['name'],x=seed['xyz'][0],y=seed['xyz'][1],z=seed['xyz'][2],
                             nfix=0,nep=0,mad_mm=0.0,ant=h['ant'],hgt=h['hgt'],seed=1,km=0.0)); continue
        pos=os.path.join(a.out,'%s.pos'%h['name'])
        cmd=[os.path.join(a.bin,'rnx2rtkp'),'-k',conf,'-o',pos]
        if a.hours and h['date']:
            y,m,d=h['date'][:3]
            cmd+=['-ts','%04d/%02d/%02d'%(y,m,d),'%02d:00:00'%a.hours[0],
                  '-te','%04d/%02d/%02d'%(y,m,d),'%02d:59:30'%(a.hours[1]-1)]
        cmd+=[fn,a.seed]+a.nav
        r=subprocess.run(cmd,capture_output=True,text=True)
        sol=read_pos(pos) if os.path.exists(pos) else []
        if a.mode=='static-iflc':
            # combined forward/backward static float: every epoch carries the
            # smoothed constant position; use the second half of the day
            fixed=[s for s in sol if s[4] in (1,2)]
            fixed=fixed[len(fixed)//2:]
        elif a.mode in ('static-fix','static-cont'):
            fixed=[s for s in sol if s[4]==1]
            fixed=fixed[len(fixed)//2:]
        else:
            fixed=[s for s in sol if s[4]==1]
        if len(fixed)<100:
            print("# %s: only %d fixed epochs of %d %s"%(h['name'],len(fixed),len(sol),r.stderr.strip()[-200:]),file=sys.stderr)
            continue
        xyz=[statistics.median(s[i] for s in fixed) for i in (1,2,3)]
        mad=[1.4826*statistics.median(abs(s[i]-xyz[i-1]) for s in fixed) for i in (1,2,3)]
        km=math.dist(xyz,seed['xyz'])/1e3
        print("# %s: %d/%d fixed, %.1f km, MAD %.1f/%.1f/%.1f mm, header diff %.2f m"%(
            h['name'],len(fixed),len(sol),km,*[1e3*m for m in mad],math.dist(xyz,h['xyz'])),file=sys.stderr)
        rows.append(dict(id=h['name'],x=xyz[0],y=xyz[1],z=xyz[2],nfix=len(fixed),nep=len(sol),
                         mad_mm=1e3*max(mad),ant=h['ant'],hgt=h['hgt'],seed=0,km=km))
    with open(os.path.join(a.out,'coords.csv'),'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['id','x','y','z','nfix','nep','mad_mm','ant','hgt','seed','km'])
        w.writeheader()
        for r in rows: w.writerow(r)
    print("wrote %s/coords.csv (%d stations)"%(a.out,len(rows)))

def merge(a):
    days={}
    for fn in a.files:
        for r in csv.DictReader(open(fn)):
            days.setdefault(r['id'],[]).append(r)
    out=[]
    for sid,rs in days.items():
        xyz=[statistics.median(float(r[k]) for r in rs) for k in 'xyz']
        spread=[max(abs(float(r[k])-xyz[i]) for r in rs) for i,k in enumerate('xyz')]
        print("# %s: %d days, day-to-day max dev %.1f/%.1f/%.1f mm, MAD %s mm, fixed %s"%(sid,len(rs),
              *[1e3*s for s in spread],'/'.join('%.0f'%float(r['mad_mm']) for r in rs),
              '/'.join('%s'%r['nfix'] for r in rs)),file=sys.stderr)
        out.append((sid,xyz,rs[0]))
    if a.yaml:
        print("stations:")
        for sid,xyz,r in out:
            print("  - id: %s\n    ecef: [%.4f, %.4f, %.4f]\n    antenna: \"%s\"\n    arp_height_m: %.4f"%(
                sid,xyz[0],xyz[1],xyz[2],r['ant'],float(r['hgt'])))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('-b','--bin',default='build/bin')
    ap.add_argument('--atx'); ap.add_argument('--nav',nargs='*',default=[]); ap.add_argument('--seed')
    ap.add_argument('--seed-xyz',nargs=3,type=float,metavar=('X','Y','Z'),
                    help='use this ECEF position for the seed instead of its header (chaining: a '
                         'station solved from the first seed becomes the seed of its neighbours)')
    ap.add_argument('--elmask',type=float,default=10.0)
    ap.add_argument('--navsys',type=int,default=1,help='RTKLIB navsys mask (1 GPS, 8 GAL)')
    ap.add_argument('--hours',nargs=2,type=int,metavar=('H0','H1'),
                    help='use only this UTC window [H0,H1) of the day, e.g. 16 23 for the quiet '
                         'night-time ionosphere in East Asia')
    ap.add_argument('--mode',choices=['kinematic-fix','static-iflc','static-fix','static-cont'],default='static-cont',
                    help='static-iflc: 24 h static ionosphere-free float (default, mm-level on 20-50 km); '
                         'kinematic-fix: median of fixed kinematic epochs')
    ap.add_argument('-o','--out',default='results/station_coords')
    ap.add_argument('--merge',action='store_true'); ap.add_argument('--yaml',action='store_true')
    ap.add_argument('files',nargs='+')
    a=ap.parse_args()
    if a.merge: merge(a)
    else:
        if not (a.atx and a.nav and a.seed): ap.error('--atx, --nav and --seed are required')
        solve(a)
if __name__=='__main__':
    main()
