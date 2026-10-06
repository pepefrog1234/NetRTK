#!/usr/bin/env python3
"""End-to-end regression: simulate a network (15 min, quiet ionosphere),
generate a leave-one-out VRS with nrtk-gen, and check observation-level and
position-level statistics against fixed thresholds.

usage: run_regression.py <build/bin> <source dir> <work dir>"""
import json, os, subprocess, sys

bindir, srcdir, work = sys.argv[1:4]
os.makedirs(work, exist_ok=True)
scen = open(os.path.join(srcdir, 'tests/regression/scenario.yaml.in')).read().replace('@OUT@', work)
open(os.path.join(work, 'scenario.yaml'), 'w').write(scen)
subprocess.run([os.path.join(bindir, 'nrtk-sim'), os.path.join(work, 'scenario.yaml')], check=True,
               capture_output=True)
conf = open(os.path.join(work, 'sim/network.yaml')).read()
conf += '  out_dir: %s/gen\n  vrs:\n    - {name: LOO_TWE5, leave_out: TWE5}\n' % work
open(os.path.join(work, 'network.yaml'), 'w').write(conf)
subprocess.run([os.path.join(bindir, 'nrtk-gen'), '-c', os.path.join(work, 'network.yaml'), '-v', '0'],
               check=True, capture_output=True)
r = subprocess.run([sys.executable, os.path.join(srcdir, 'tools/loo_validate.py'), '-c',
                    os.path.join(work, 'network.yaml'), '-j', 'LOO_TWE5', '-b', bindir, '--reset', '300',
                    '--warmup', '120'], capture_output=True, text=True, cwd=srcdir)
print(r.stdout[-3000:])
s = json.load(open(os.path.join(work, 'gen/validation_LOO_TWE5/summary.json')))
v = s['results']['VRS']
l1 = s['obs_bands']['L1/E1/B1']
checks = {
    'network fix rate (all baselines) >= 0.85': all(float(b['fix_rate']) >= 0.85 for b in
        __import__('csv').DictReader(open(os.path.join(work, 'gen/baselines.csv')))),
    'L1 DD phase residual <= 0.015 m': l1['phase_rms'] <= 0.015,
    'VRS rover fix rate >= 0.95': v['fix_rate'] >= 0.95,
    'VRS horizontal RMS <= 0.02 m': v.get('h_rms', 1) <= 0.02,
    'VRS vertical RMS <= 0.04 m': v.get('v_rms', 1) <= 0.04,
    'no wrong fixes > 10 cm': v.get('wrong_fix', 1) <= 0.001,
}
ok = True
for k, c in checks.items():
    print('[%s] %s' % ('PASS' if c else 'FAIL', k))
    ok &= bool(c)
sys.exit(0 if ok else 1)
