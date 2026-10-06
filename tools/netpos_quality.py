#!/usr/bin/env python3
"""Quality of fixed network residuals without truth (real data):
|D - DD GMF model| distribution (wide-lane errors show up at ~0.32-0.38 m,
L1 errors at ~0.11 m), ionosphere magnitudes and rates.
usage: netpos_quality.py blsat_<A>-<B>.csv [...]"""
import csv, collections, datetime, sys
def t2s(ts):
    t = datetime.datetime.strptime(ts[:19], '%Y/%m/%d %H:%M:%S'); return t.hour*3600 + t.minute*60 + t.second
for fn in sys.argv[1:]:
    rows = list(csv.DictReader(open(fn)))
    refs = {(r['time'], r['sat'][0]): r['sat'] for r in rows if r['isref'] == '1'}
    nd, io, n = [], [], 0
    ser = collections.defaultdict(list)
    for r in rows:
        if r['fixed'] != '1' or r['isref'] == '1':
            continue
        n += 1
        nd.append(abs(float(r['dd_nondisp']) - float(r['dd_trop_model'])))
        io.append(abs(float(r['dd_iono'])))
    if not nd:
        print(fn, 'no fixed'); continue
    nd.sort(); io.sort()
    q = lambda v, p: v[min(len(v)-1, int(p*len(v)))]
    print('%-34s n=%6d |D-model| p50 %.3f p95 %.3f p99 %.3f >0.15m %.2f%% | |I| p50 %.3f p95 %.3f max %.3f' % (
        fn.split('/')[-1], n, q(nd, .5), q(nd, .95), q(nd, .99), 100.0*sum(1 for x in nd if x > 0.15)/len(nd), q(io, .5), q(io, .95), io[-1]))
