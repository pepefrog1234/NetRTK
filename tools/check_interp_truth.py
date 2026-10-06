#!/usr/bin/env python3
"""Compare interpolated VRS corrections (nrtk-gen *_diag.csv) with simulator truth.

usage: check_interp_truth.py job_diag.csv truth_dir VRS_STATION
Truth DD (master -> VRS station) = (X_V^s - X_M^s) - (X_V^r - X_M^r).
"""
import csv, sys, math, collections, datetime

GPS0 = datetime.datetime(1980, 1, 6)

def tow_of(ts):
    t = datetime.datetime.strptime(ts[:19], '%Y/%m/%d %H:%M:%S')
    return round(((t - GPS0).total_seconds() + float('0' + ts[19:])) % 604800.0, 1)

def load_truth(fn):
    d = {}
    for r in csv.DictReader(open(fn)):
        d[(float(r['tow']), r['sat'])] = (float(r['iono_l1']), float(r['trop']) + float(r['orbit']), float(r['el']))
    return d

diag, tdir, vsta = sys.argv[1:4]
TV = load_truth('%s/truth_%s.csv' % (tdir, vsta))
TM = {}
rows = collections.defaultdict(list)
for r in csv.DictReader(open(diag)):
    rows[tow_of(r['time'])].append(r)
stats = collections.defaultdict(lambda: [[], []])
worst = []
for tow, rs in rows.items():
    if (tow, rs[0]['sat']) not in TV and not any((tow, r['sat']) in TV for r in rs):
        continue
    refs = {r['sat'][0]: r['sat'] for r in rs if r['isref'] == '1'}
    for r in rs:
        if r['ok'] != '1' or r['isref'] == '1':
            continue
        m = r['master']
        if m not in TM:
            TM[m] = load_truth('%s/truth_%s.csv' % (tdir, m))
        s, ref = r['sat'], refs[r['sat'][0]]
        k1, k2 = (tow, s), (tow, ref)
        if not (k1 in TV and k2 in TV and k1 in TM[m] and k2 in TM[m]):
            continue
        ti = (TV[k1][0] - TM[m][k1][0]) - (TV[k2][0] - TM[m][k2][0])
        tn = (TV[k1][1] - TM[m][k1][1]) - (TV[k2][1] - TM[m][k2][1])
        ei, en = float(r['dd_iono']) - ti, float(r['dd_nondisp']) - tn
        el = TV[k1][2]
        key = 'nbl=%s extrap=%s el%s' % (r['nbl'], r['extrap'], '<20' if el < 20 else ('20-40' if el < 40 else '>=40'))
        stats[key][0].append(ei); stats[key][1].append(en)
        worst.append((abs(ei - 0) + abs(en), r['time'], s, ref, r['nbl'], r['extrap'], round(el, 1), round(ei, 4), round(en, 4)))
rms = lambda v: math.sqrt(sum(x * x for x in v) / len(v)) if v else 0.0
print('%-28s %6s %9s %9s %9s %9s %10s' % ('group', 'n', 'I_mean', 'I_rms', 'D_mean', 'D_rms', 'L1err_rms'))
for k in sorted(stats):
    a, b = stats[k]
    l1 = [y - x for x, y in zip(a, b)]  # L1 phase correction error = D - I
    print('%-28s %6d %9.4f %9.4f %9.4f %9.4f %10.4f' % (k, len(a), sum(a)/len(a), rms(a), sum(b)/len(b), rms(b), rms(l1)))
worst.sort(reverse=True)
print('worst:', worst[:5])
