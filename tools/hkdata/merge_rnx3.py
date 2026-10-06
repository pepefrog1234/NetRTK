#!/usr/bin/env python3
"""Merge consecutive RINEX 3 observation files of ONE station (e.g. hourly -> daily).
- header taken from the first file; '# OF SATELLITES' / 'PRN / # OF OBS' lines dropped; TIME OF LAST OBS updated
- refuses to merge if SYS / # / OBS TYPES differ between inputs
- drops epoch-flag-4 'header information' event blocks (receiver statistics) found inside the inputs
usage: merge_rnx3.py out.rnx in1.rnx in2.rnx ..."""
import sys
def split(fn):
    hdr=[]; body=[]
    with open(fn) as f:
        for line in f:
            hdr.append(line)
            if line[60:].startswith('END OF HEADER'): break
        body=f.readlines()
    return hdr, body
def obstypes(hdr): return [l for l in hdr if l[60:].startswith('SYS / # / OBS TYPES')]
def clean_body(body):
    out=[]; i=0; n=len(body)
    while i<n:
        l=body[i]
        if l.startswith('>'):
            flag=int(l[31]); ns=int(l[32:35])
            if flag in (4,):  # header-info event block -> drop it and its ns lines
                i+=1+ns; continue
            out.append(l); out.extend(body[i+1:i+1+ns]); i+=1+ns; continue
        i+=1  # stray line (should not happen)
    return out
def main():
    out=sys.argv[1]; ins=sys.argv[2:]
    h0,b0=split(ins[0]); ot0=obstypes(h0)
    bodies=[clean_body(b0)]
    for fn in ins[1:]:
        h,b=split(fn)
        if obstypes(h)!=ot0:
            sys.exit('OBS TYPES differ in %s - not merging'%fn)
        bodies.append(clean_body(b))
    last=None
    for b in reversed(bodies):
        for l in reversed(b):
            if l.startswith('>') and l[31]=='0' or (l.startswith('>') and l[31]=='1'):
                last=l; break
        if last: break
    hdr=[]
    for l in h0:
        lab=l[60:].rstrip()
        if lab in ('# OF SATELLITES','PRN / # OF OBS'): continue
        if lab=='TIME OF LAST OBS' and last:
            y,mo,d,hh,mi=[int(x) for x in last[2:29].split()[:5]]; s=float(last[2:29].split()[5])
            l='%6d%6d%6d%6d%6d%13.7f%8s%9s%s\n'%(y,mo,d,hh,mi,s,'GPS','','TIME OF LAST OBS')
            l=l[:60].ljust(60)+'TIME OF LAST OBS\n'
        hdr.append(l)
    # if header had no TIME OF LAST OBS, leave as is
    with open(out,'w') as f:
        f.writelines(hdr)
        for b in bodies: f.writelines(b)
    print('merged %d files -> %s (last epoch %s)'%(len(ins),out,last[2:29].strip() if last else '?'))
if __name__=='__main__': main()
