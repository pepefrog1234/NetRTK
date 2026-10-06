#!/usr/bin/env python3
"""Does the LSM fit residual (rms_iono) predict the true interpolation error?
usage: check_quality_indicator.py job_diag.csv truth_dir VRS_STATION"""
import csv, collections, datetime, sys
GPS0 = datetime.datetime(1980, 1, 6)
def tow_of(ts):
    t = datetime.datetime.strptime(ts[:19], '%Y/%m/%d %H:%M:%S'); return round(((t - GPS0).total_seconds()) % 604800.0, 1)
def load(fn):
    return {(float(r['tow']), r['sat']): float(r['iono_l1']) for r in csv.DictReader(open(fn))}
diag, tdir, vsta = sys.argv[1:4]
TV, TM = load('%s/truth_%s.csv' % (tdir, vsta)), {}
rows = collections.defaultdict(list)
for r in csv.DictReader(open(diag)):
    if r['time'][17:19] in ('00', '10', '20', '30', '40', '50'):
        rows[tow_of(r['time'])].append(r)
pairs = []
for tow, rs in rows.items():
    refs = {r['sat'][0]: r['sat'] for r in rs if r['isref'] == '1'}
    for r in rs:
        if r['ok'] != '1' or r['isref'] == '1' or r['nbl'] != '3':
            continue
        m = r['master']
        if m not in TM:
            TM[m] = load('%s/truth_%s.csv' % (tdir, m))
        k1, k2 = (tow, r['sat']), (tow, refs[r['sat'][0]])
        if not all(k in TV and k in TM[m] for k in (k1, k2)):
            continue
        ti = (TV[k1] - TM[m][k1]) - (TV[k2] - TM[m][k2])
        pairs.append((float(r['rms_iono']), abs(float(r['dd_iono']) - ti)))
pairs.sort(); n = len(pairs)
for q in (0.5, 0.75, 0.9, 0.95):
    thr = pairs[int(q * n) - 1][0]
    hi = [e for x, e in pairs if x > thr]; lo = [e for x, e in pairs if x <= thr]
    print('rms_iono > %.4f m (top %2.0f%%): mean |err| %.4f m vs %.4f m below' % (thr, 100 * (1 - q), sum(hi) / len(hi), sum(lo) / len(lo)))
xs = [p[0] for p in pairs]; ys = [p[1] for p in pairs]; mx, my = sum(xs) / n, sum(ys) / n
c = sum((x - mx) * (y - my) for x, y in pairs) / n
sx = (sum((x - mx) ** 2 for x in xs) / n) ** .5; sy = (sum((y - my) ** 2 for y in ys) / n) ** .5
print('corr(rms_iono, |iono error|) = %.2f  (n=%d)' % (c / sx / sy, n))
