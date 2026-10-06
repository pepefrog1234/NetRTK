#!/usr/bin/env python3
"""Leave-one-out validation (spec ch.10 / ch.12).

The left-out station is processed as a rover with RTKLIB rnx2rtkp against
  (1) the VRS generated at its coordinates (nrtk-gen job),
  (2) the nearest physical station of the remaining network,
  (3) the master station mostly used by the VRS,
and the solutions are compared with its known coordinates. The observation
layer (VRS vs station double differences) is taken from `nrtk-tools obscmp`,
the network layer from nrtk-gen's baselines.csv.

usage:
  loo_validate.py -c config.yaml -j JOB [-b build/bin] [-o report_dir]
                  [--reset 600] [--systems GEJ] [--elmask 10]
"""
import argparse, csv, datetime, json, math, os, subprocess, sys, collections

GPS0 = datetime.datetime(1980, 1, 6)
SYSMASK = {'G': 1, 'R': 4, 'E': 8, 'J': 16, 'C': 32}


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.stderr.write(' '.join(cmd) + '\n' + r.stdout[-2000:] + r.stderr[-2000:])
        raise SystemExit('command failed')
    return r.stdout


def ecef2llh(x, y, z):
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    r2 = x * x + y * y
    zk, v, z0 = 0.0, a, z
    while True:
        zk = z0
        sinp = z0 / math.sqrt(r2 + z0 * z0)
        v = a / math.sqrt(1 - e2 * sinp * sinp)
        z0 = z + v * e2 * sinp
        if abs(z0 - zk) < 1e-4:
            break
    lat = math.atan(z0 / math.sqrt(r2)) if r2 > 1e-12 else (math.pi / 2 if z > 0 else -math.pi / 2)
    lon = math.atan2(y, x) if r2 > 1e-12 else 0.0
    return lat, lon, math.sqrt(r2 + z0 * z0) - v


def enu(ref, p):
    lat, lon, _ = ecef2llh(*ref)
    d = [p[i] - ref[i] for i in range(3)]
    sl, cl, sp, cp = math.sin(lon), math.cos(lon), math.sin(lat), math.cos(lat)
    e = -sl * d[0] + cl * d[1]
    n = -sp * cl * d[0] - sp * sl * d[1] + cp * d[2]
    u = cp * cl * d[0] + cp * sl * d[1] + sp * d[2]
    return e, n, u


def rtk_conf(path, antex, rover, base, base_pos, base_ant, base_del, navsys, elmask, iono):
    lines = {
        'pos1-posmode': 'kinematic', 'pos1-frequency': 'l1+l2', 'pos1-soltype': 'forward',
        'pos1-elmask': str(elmask), 'pos1-dynamics': 'off', 'pos1-tidecorr': 'off',
        'pos1-ionoopt': iono, 'pos1-tropopt': 'saas', 'pos1-sateph': 'brdc',
        'pos1-navsys': str(navsys), 'pos2-armode': 'continuous', 'pos2-gloarmode': 'off',
        'pos2-bdsarmode': 'on', 'pos2-arthres': '3', 'pos2-arfilter': 'on',
        'pos2-elmaskar': '15', 'pos2-maxage': '30',
        'out-solformat': 'xyz', 'out-outhead': 'on', 'out-timeform': 'tow', 'out-timesys': 'gpst',
        'ant1-postype': 'single', 'ant1-anttype': rover['antenna'] or '',
        'ant1-antdele': str(rover['antdel'][0]), 'ant1-antdeln': str(rover['antdel'][1]),
        'ant1-antdelu': str(rover['antdel'][2]),
        'ant2-postype': 'xyz', 'ant2-pos1': '%.4f' % base_pos[0], 'ant2-pos2': '%.4f' % base_pos[1],
        'ant2-pos3': '%.4f' % base_pos[2], 'ant2-anttype': base_ant or '',
        'ant2-antdele': str(base_del[0]), 'ant2-antdeln': str(base_del[1]), 'ant2-antdelu': str(base_del[2]),
        'file-rcvantfile': antex,
    }
    with open(path, 'w') as f:
        for k, v in lines.items():
            f.write('%-20s=%s\n' % (k, v))


def read_pos(fn):
    out = []
    for line in open(fn):
        if line.startswith('%') or not line.strip():
            continue
        p = line.split()
        if len(p) < 6:
            continue
        week, tow = int(p[0]), float(p[1])
        out.append((week * 604800.0 + tow, float(p[2]), float(p[3]), float(p[4]), int(p[5]), int(p[6])))
    return out


def stats(sol, truth, t0=None):
    n = len(sol)
    fix = [s for s in sol if s[4] == 1]
    res = {'epochs': n, 'fixed': len(fix), 'fix_rate': len(fix) / n if n else 0.0}
    if fix:
        d = [enu(truth, s[1:4]) for s in fix]
        h = [math.hypot(e, nn) for e, nn, u in d]
        v = [u for e, nn, u in d]
        d3 = [math.sqrt(e * e + nn * nn + u * u) for e, nn, u in d]
        res['h_rms'] = math.sqrt(sum(x * x for x in h) / len(h))
        res['v_rms'] = math.sqrt(sum(x * x for x in v) / len(v))
        res['v_mean'] = sum(v) / len(v)
        res['max3d'] = max(d3)
        res['wrong_fix'] = sum(1 for x in d3 if x > 0.10) / len(d3)
    return res


def hourly(sol, truth):
    """per UTC hour statistics of a solution"""
    out = {}
    for s in sol:
        h = int((s[0] % 86400) // 3600)
        out.setdefault(h, []).append(s)
    return {h: stats(v, truth) for h, v in sorted(out.items())}


def gps_seconds(s):
    t = datetime.datetime.strptime(s, '%Y/%m/%d %H:%M:%S')
    return (t - GPS0).total_seconds()


def solve(bin_dir, conf, rover_file, base_file, nav, out, ts=None, te=None):
    cmd = [os.path.join(bin_dir, 'rnx2rtkp'), '-k', conf, '-o', out]
    if ts:
        cmd += ['-ts'] + ts.split()
    if te:
        cmd += ['-te'] + te.split()
    cmd += [rover_file, base_file] + nav.split()
    run(cmd)
    return read_pos(out)


def ttff(bin_dir, conf, rover_file, base_file, nav, outdir, tag, start, end, reset):
    t, res = gps_seconds(start), []
    tend = gps_seconds(end)
    k = 0
    while t + reset <= tend + 1e-6:
        ts = (GPS0 + datetime.timedelta(seconds=t)).strftime('%Y/%m/%d %H:%M:%S')
        te = (GPS0 + datetime.timedelta(seconds=t + reset - 1)).strftime('%Y/%m/%d %H:%M:%S')
        sol = solve(bin_dir, conf, rover_file, base_file, nav, os.path.join(outdir, '%s_seg%02d.pos' % (tag, k)), ts, te)
        first = next((s[0] for s in sol if s[4] == 1), None)
        res.append(first - t if first is not None else float('nan'))
        t += reset
        k += 1
    return res


def obs_layer(bin_dir, conf, sta, vrs_obs, out_csv):
    txt = run([os.path.join(bin_dir, 'nrtk-tools'), 'obscmp', '-c', conf, '-s', sta, '-v', vrs_obs, '-o', out_csv])
    bands = {}
    section = None
    for line in txt.splitlines():
        if line.startswith('band:'):
            section = 'band'
            continue
        p = line.split()
        if section == 'band' and len(p) == 5:
            bands[p[0]] = {'n': int(p[1]), 'phase_rms': float(p[2]), 'phase_max': float(p[3]), 'code_rms': float(p[4])}
    return txt, bands


def plot(outdir, sols, truth, title):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    fig, axs = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    colors = {'VRS': '#1f77b4', 'nearest': '#d62728', 'master': '#7f7f7f'}
    for name, sol in sols.items():
        if not sol:
            continue
        t0 = sol[0][0]
        for k, lab in enumerate(['E', 'N', 'U']):
            xs = [(s[0] - t0) / 60.0 for s in sol if s[4] == 1]
            ys = [enu(truth, s[1:4])[k] * 100 for s in sol if s[4] == 1]
            axs[k].plot(xs, ys, '.', ms=1.5, color=colors.get(name, None), label=name if k == 0 else None)
            axs[k].set_ylabel('d%s (cm)' % lab)
            axs[k].set_ylim(-10, 10)
            axs[k].grid(alpha=0.3)
    axs[0].legend(markerscale=6, loc='upper right')
    axs[2].set_xlabel('minutes (fixed epochs only)')
    fig.suptitle(title)
    fn = os.path.join(outdir, 'enu_fixed.png')
    fig.tight_layout()
    fig.savefig(fn, dpi=120)
    return fn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-c', '--conf', required=True)
    ap.add_argument('-j', '--job', required=True)
    ap.add_argument('-b', '--bin', default='build/bin')
    ap.add_argument('-o', '--out', default=None)
    ap.add_argument('--reset', type=float, default=600.0, help='TTFF reset interval (s)')
    ap.add_argument('--systems', default='GEJC')
    ap.add_argument('--elmask', type=float, default=10.0)
    ap.add_argument('--iono', default='off', help='rover ionosphere option for all three solutions')
    ap.add_argument('--warmup', type=float, default=120.0,
                    help='skip the network cold start (s) for all rover solutions')
    a = ap.parse_args()

    cfg = json.loads(run([os.path.join(a.bin, 'nrtk-tools'), 'confjson', a.conf]))
    job = next(j for j in cfg['jobs'] if j['name'] == a.job)
    sta = {s['id']: s for s in cfg['stations']}
    loo = sta[job['leave_out']]
    gen_dir = cfg['outdir']
    out = a.out or os.path.join(gen_dir, 'validation_' + a.job)
    os.makedirs(out, exist_ok=True)
    vrs_obs = os.path.join(gen_dir, a.job + '.obs')
    navsys = sum(SYSMASK[c] for c in a.systems)
    truth = loo['ecef']
    vrs_pos = truth if job['leave_out'] else job['ecef']

    # master most used by the VRS; nearest physical station (excluding the left-out one)
    masters = collections.Counter(r['master'] for r in csv.DictReader(open(os.path.join(gen_dir, a.job + '_epoch.csv'))) if r['mode'] == '1')
    master = masters.most_common(1)[0][0]
    others = [s for s in cfg['stations'] if s['id'] != loo['id'] and s['enabled']]
    nearest = min(others, key=lambda s: math.dist(s['ecef'], truth))

    cases = {
        'VRS': (vrs_obs, vrs_pos, '', [0.0, 0.0, 0.0]),
        'nearest': (nearest['rinex'], nearest['ecef'], nearest['antenna'], nearest['antdel']),
        'master': (sta[master]['rinex'], sta[master]['ecef'], sta[master]['antenna'], sta[master]['antdel']),
    }
    results, sols, ttffs = {}, {}, {}
    t_start = None
    if cfg['start'][:4] != '1970':
        t_start = (GPS0 + datetime.timedelta(seconds=gps_seconds(cfg['start']) + a.warmup)).strftime('%Y/%m/%d %H:%M:%S')
    t_end = cfg['end'] if cfg['end'][:4] != '1970' else None
    for name, (bfile, bpos, bant, bdel) in cases.items():
        if name == 'master' and master == nearest['id']:
            continue
        conf = os.path.join(out, 'rtk_%s.conf' % name)
        rtk_conf(conf, cfg['antex'], loo, bfile, bpos, bant, bdel, navsys, a.elmask, a.iono)
        sol = solve(a.bin, conf, loo['rinex'], bfile, cfg['nav'], os.path.join(out, '%s.pos' % name),
                    t_start, t_end)
        sols[name] = sol
        results[name] = stats(sol, truth)
        results[name]['hourly'] = hourly(sol, truth)
        results[name]['base_km'] = math.dist(bpos, truth) / 1e3
        if a.reset > 0 and sol:
            t0 = (GPS0 + datetime.timedelta(seconds=sol[0][0])).strftime('%Y/%m/%d %H:%M:%S')
            t1 = (GPS0 + datetime.timedelta(seconds=sol[-1][0])).strftime('%Y/%m/%d %H:%M:%S')
            tt = ttff(a.bin, conf, loo['rinex'], bfile, cfg['nav'], out, name, t0, t1, a.reset)
            ok = [x for x in tt if x == x]
            ttffs[name] = tt
            results[name]['ttff_mean'] = sum(ok) / len(ok) if ok else float('nan')
            results[name]['ttff_fail'] = len(tt) - len(ok)

    obs_txt, bands = obs_layer(a.bin, a.conf, loo['id'], vrs_obs, os.path.join(out, 'obscmp.csv'))
    bl = list(csv.DictReader(open(os.path.join(gen_dir, 'baselines.csv')))) if os.path.exists(os.path.join(gen_dir, 'baselines.csv')) else []

    # acceptance (spec ch.12, prototype) -------------------------------------
    v, nr = results.get('VRS', {}), results.get('nearest', {})
    l1 = bands.get('L1/E1/B1', {})
    checks = [
        ('Network: baseline fix rate >= 95% (all baselines)', all(float(b['fix_rate']) >= 0.95 for b in bl) if bl else None,
         ', '.join('%s %.1f%%' % (b['baseline'], 100 * float(b['fix_rate'])) for b in bl)),
        ('Network: cold-start TTFF <= 120 s', all(0 <= float(b['ttff_s']) <= 120 for b in bl) if bl else None,
         ', '.join('%s %ss' % (b['baseline'], b['ttff_s']) for b in bl)),
        ('Observation: L1 DD phase residual RMS <= 2 cm (el >= 15)', l1.get('phase_rms', 9) <= 0.02, '%.2f cm' % (100 * l1.get('phase_rms', float('nan')))),
        ('Observation: code residual RMS <= 0.5 m', l1.get('code_rms', 9) <= 0.5, '%.2f m' % l1.get('code_rms', float('nan'))),
        ('Position: VRS horizontal RMS <= 2 cm (fixed)', v.get('h_rms', 9) <= 0.02, '%.2f cm' % (100 * v.get('h_rms', float('nan')))),
        ('Position: VRS vertical RMS <= 4 cm (fixed)', v.get('v_rms', 9) <= 0.04, '%.2f cm' % (100 * v.get('v_rms', float('nan')))),
        ('Position: VRS fix rate >= 95% and >= nearest + 5 pts', v.get('fix_rate', 0) >= 0.95 and v.get('fix_rate', 0) >= nr.get('fix_rate', 0) + 0.05,
         'VRS %.1f%%, nearest %.1f%%' % (100 * v.get('fix_rate', 0), 100 * nr.get('fix_rate', 0))),
        ('Position: rover TTFF <= 60 s and < nearest', v.get('ttff_mean', 1e9) <= 60 and v.get('ttff_mean', 1e9) < nr.get('ttff_mean', 1e9) if nr.get('ttff_mean', 1e9) == nr.get('ttff_mean', 1e9) else v.get('ttff_mean', 1e9) <= 60,
         'VRS %.0f s, nearest %s' % (v.get('ttff_mean', float('nan')),
                                     '%.0f s' % nr['ttff_mean'] if nr.get('ttff_mean') is not None and nr['ttff_mean'] == nr['ttff_mean'] else 'no fix')),
        ('Position: wrong fixes (>10 cm) <= 0.1%', v.get('wrong_fix', 1) <= 0.001, '%.3f%%' % (100 * v.get('wrong_fix', float('nan')))),
    ]
    png = plot(out, sols, truth, '%s: leave-one-out %s' % (cfg['name'], loo['id']))
    summary = {'job': a.job, 'station': loo['id'], 'master': master, 'nearest': nearest['id'],
               'results': results, 'ttff_segments': ttffs, 'obs_bands': bands,
               'checks': [{'item': c[0], 'pass': c[1], 'value': c[2]} for c in checks]}
    json.dump(summary, open(os.path.join(out, 'summary.json'), 'w'), indent=1, default=str)

    print('Leave-one-out %s (%s): master %s, nearest %s' % (loo['id'], a.job, master, nearest['id']))
    print('%-8s %7s %8s %8s %8s %8s %8s %8s %6s' % ('base', 'km', 'fix%', 'H_rms', 'V_rms', 'max3d', 'wrong%', 'TTFF', 'fail'))
    for name, r in results.items():
        print('%-8s %7.1f %7.1f%% %7.2fcm %7.2fcm %7.2fcm %7.3f%% %7.0fs %6s' % (
            name, r['base_km'], 100 * r['fix_rate'], 100 * r.get('h_rms', float('nan')), 100 * r.get('v_rms', float('nan')),
            100 * r.get('max3d', float('nan')), 100 * r.get('wrong_fix', float('nan')), r.get('ttff_mean', float('nan')),
            r.get('ttff_fail', '')))
    print(obs_txt)
    if len(results.get('VRS', {}).get('hourly', {})) > 2:
        print('hour(UTC)  ' + '  '.join('%-24s' % n for n in results))
        hours = sorted(set(h for r in results.values() for h in r['hourly']))
        for h in hours:
            row = []
            for n, r in results.items():
                x = r['hourly'].get(h)
                row.append('fix %5.1f%% H %5.1fcm w%5.2f%%' % (100 * x['fix_rate'], 100 * x.get('h_rms', float('nan')),
                                                         100 * x.get('wrong_fix', float('nan'))) if x else '-')
            print('%02d         %s' % (h, '  '.join('%-24s' % c for c in row)))
    for item, ok, val in checks:
        print('[%s] %s: %s' % ('PASS' if ok else ('n/a' if ok is None else 'FAIL'), item, val))
    if png:
        print('plot:', png)


if __name__ == '__main__':
    main()
