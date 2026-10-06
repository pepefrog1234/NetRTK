#!/usr/bin/env python3
"""Derive the configuration of another day from a CWA/GDMS day configuration.

  tools/twdata/config_for_day.py conf/tw/tw_north_30s.yaml 2026 271 > conf/tw/tw_north_30s_271.yaml

Rewrites the RINEX paths (data/real/tw/YYYY_DDD/ssssDDD0.YYo.gz), the broadcast
navigation file, offline.start and out_dir. Station coordinates stay the same:
processing a second day with them is the independent check of their quality
(tools/twdata/adjust_coords.py should then report corrections at the mm level).
"""
import datetime, re, sys

def main():
    if len(sys.argv)!=4: sys.exit(__doc__)
    conf,year,doy=sys.argv[1],int(sys.argv[2]),int(sys.argv[3])
    txt=open(conf).read()
    yy='%02d'%(year%100)
    date=datetime.date(year,1,1)+datetime.timedelta(days=doy-1)
    txt=re.sub(r'data/real/tw/%d_\d{3}/([a-z0-9]{4})\d{3}0\.%so\.gz'%(year,yy),
               r'data/real/tw/%d_%03d/\g<1>%03d0.%so.gz'%(year,doy,doy,yy),txt)
    txt=re.sub(r'BRDC00IGS_R_%d\d{3}0000_01D_MN\.rnx'%year,'BRDC00IGS_R_%d%03d0000_01D_MN.rnx'%(year,doy),txt)
    txt=re.sub(r'start: "\d{4}-\d{2}-\d{2} 00:00:00"','start: "%s 00:00:00"'%date.isoformat(),txt)
    txt=re.sub(r'(out_dir: \S+?)(_\d{3})?\n',r'\g<1>_%03d\n'%doy,txt,count=1)
    sys.stdout.write(txt)
if __name__=='__main__':
    main()
