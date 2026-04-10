#!/usr/bin/env python3
"""Quick accuracy test for all inelastic reactions.
Run from Cpp_AI directory: python3 test_all_inel.py
"""
import subprocess, numpy as np, re, os

CPPDIR = os.path.dirname(os.path.abspath(__file__))

def parse_cleopatra(fn):
    r = {}; in_dcs = False
    pat = re.compile(r'^\s{1,6}(\d{1,3}\.\d{2})\s+(0\.\d+E[+-]\d{2}|\d+\.?\d*)')
    for line in open(fn):
        if 'COMPUTATION OF CROSS SECTIONS' in line: in_dcs = True
        if 'ANALYZING POWERS' in line: in_dcs = False
        if not in_dcs: continue
        m = pat.match(line)
        if not m: continue
        try:
            a = float(m.group(1)); d = float(m.group(2).rstrip('.'))
            if 0 <= a <= 180 and 0 < d < 1e15:
                r[round(a, 1)] = d
        except: pass
    return r

def run_cpp(inp):
    r = subprocess.run(['./ptolemy_cpp'], stdin=open(inp),
                       capture_output=True, text=True, timeout=30, cwd=CPPDIR)
    cd = {}
    for line in r.stdout.split('\n'):
        parts = line.split()
        if len(parts) >= 2:
            try:
                a = float(parts[0]); d = float(parts[1])
                if d > 0: cd[round(a, 1)] = d
            except: pass
    return cd

def err(cpp, ref):
    a = sorted(set(cpp) & set(ref))
    e = [abs(cpp[v]-ref[v])/ref[v]*100 for v in a if ref.get(v, 0) > 1e-30]
    if not e: return 0, 0, 0
    se = sorted(enumerate(e), key=lambda x: -x[1])
    return np.mean(e), se[0][1], a[se[0][0]]

tests = [
    ('206Hg(d,d\') 4+', 'test_hg206dd_inel.in', 'fortran_inel_ref.dat', False),
    ('16O(p,p\') 3-',   'test_inputs/o16_pp_3minus.in', 'test_inputs/ftn_o16pp.txt', True),
    ('32Si(d,d\') 2+',  'test_inputs/si32_dd_inel.in',  'test_inputs/ftn_si32.txt',  True),
    ('48Ca(a,a\') 2+',  'test_inputs/ca48_aa_inel.in',  'test_inputs/ftn_ca48.txt',  True),
]

print("=== Ptolemy C++ Inelastic Accuracy Test ===\n")
for label, inp, rfile, cleo in tests:
    inp_path = os.path.join(CPPDIR, inp)
    ref_path = os.path.join(CPPDIR, rfile)
    cpp = run_cpp(inp_path)
    ref = parse_cleopatra(ref_path) if cleo else {r[0]: r[1] for r in np.loadtxt(ref_path)}
    m, mx, a = err(cpp, ref)
    status = '✅' if m < 1.0 else '⚠️' if m < 10 else '❌'
    print(f"  {status} {label:20s}: Mean={m:.4f}%  Max={mx:.4f}% at {a:.0f}°")

print("\nDone.")
