#!/usr/bin/env python3
"""Compare netpos DD ionosphere / non-dispersive estimates with simulator truth.

usage: check_netpos_truth.py blsat_<M>-<A>.csv truth_<M>.csv truth_<A>.csv
DD convention: (A - M) for satellite s minus the same for the reference r.
"""
import csv, sys, math, collections, datetime

GPS0 = datetime.datetime(1980, 1, 6)

def tow_of(ts):
    t = datetime.datetime.strptime(ts[:19], '%Y/%m/%d %H:%M:%S')
    sec = (t - GPS0).total_seconds() + float('0' + ts[19:])
    return round(sec % 604800.0, 1)

def load_truth(fn):
    d = {}
    for r in csv.DictReader(open(fn)):
        d[(float(r['tow']), r['sat'])] = (float(r['iono_l1']), float(r['trop']) + float(r['orbit']))
    return d

def main():
    bl, tm, ta = sys.argv[1:4]
    TM, TA = load_truth(tm), load_truth(ta)
    rows = collections.defaultdict(list)
    for r in csv.DictReader(open(bl)):
        rows[tow_of(r['time'])].append(r)
    di, dd = collections.defaultdict(list), collections.defaultdict(list)
    for tow, rs in rows.items():
        refs = {r['sat'][0]: r['sat'] for r in rs if r['isref'] == '1'}
        for r in rs:
            if r['fixed'] != '1' or r['isref'] == '1':
                continue
            s, ref = r['sat'], refs.get(r['sat'][0])
            if ref is None:
                continue
            k = [(tow, s), (tow, ref)]
            if not all(x in TM and x in TA for x in k):
                continue
            ti = (TA[k[0]][0] - TM[k[0]][0]) - (TA[k[1]][0] - TM[k[1]][0])
            tn = (TA[k[0]][1] - TM[k[0]][1]) - (TA[k[1]][1] - TM[k[1]][1])
            el = float(r['el_m'])
            band = '10-20' if el < 20 else ('20-40' if el < 40 else '40-90')
            di[band].append(float(r['dd_iono']) - ti)
            dd[band].append(float(r['dd_nondisp']) - tn)
    print('%-7s %7s %10s %10s %10s %10s' % ('el', 'n', 'iono_mean', 'iono_rms', 'nd_mean', 'nd_rms'))
    allI, allD = [], []
    for band in sorted(di):
        a, b = di[band], dd[band]
        allI += a; allD += b
        rms = lambda v: math.sqrt(sum(x * x for x in v) / len(v))
        print('%-7s %7d %10.4f %10.4f %10.4f %10.4f' % (band, len(a), sum(a)/len(a), rms(a), sum(b)/len(b), rms(b)))
    if allI:
        bad = sum(1 for x in allD if abs(x) > 0.05)
        print('total %d, |nd err|>5cm: %d (%.2f%%), max |iono err| %.3f m, max |nd err| %.3f m' %
              (len(allI), bad, 100.0*bad/len(allI), max(abs(x) for x in allI), max(abs(x) for x in allD)))

if __name__ == '__main__':
    main()
