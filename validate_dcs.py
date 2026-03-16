#!/usr/bin/env python3
"""
validate_dcs.py — Independent DCS validation for 148Sm(α,α) at 50 MeV
Compares C++ DCS formula against Python reference implementation.
"""
from mpmath import mp, loggamma
import numpy as np

mp.dps = 25

# ── Physics parameters (must match Ptolemy / Raphael) ─────────────────────────
Ap, At, Zp, Zt = 4.0, 148.0, 2.0, 62.0
Elab = 50.0
mu   = Ap * At / (Ap + At)           # reduced mass [AMU]
Ecm  = Elab * At / (Ap + At)         # CM energy [MeV]
AMU  = 931.494061                     # MeV/c^2
HBARC = 197.32697                     # MeV·fm

k    = float(mp.sqrt(2 * mu * AMU * Ecm)) / HBARC   # fm^-1
eta  = Zp * Zt * mu * AMU / (137.036 * HBARC * k)

print(f"Kinematics check:")
print(f"  Ecm   = {Ecm:.4f} MeV  (Ptolemy: 48.683 MeV)")
print(f"  k     = {k:.4f} fm^-1 (Ptolemy: 3.0124 fm^-1)")
print(f"  eta   = {eta:.4f}      (Ptolemy: should be ~21)")
print()

# ── Coulomb phase shift σ_L = Im[ln Γ(L+1+iη)] ───────────────────────────────
def sigma_L(L):
    """Coulomb phase shift via mpmath.loggamma (exact)."""
    return float(loggamma(complex(L + 1, eta)).imag)

# ── Legendre polynomial P_L(x) via recurrence ─────────────────────────────────
def legendre(L, x):
    if L == 0: return 1.0
    if L == 1: return float(x)
    pm2, pm1 = 1.0, float(x)
    for l in range(2, L + 1):
        p = ((2*l - 1) * x * pm1 - (l - 1) * pm2) / l
        pm2, pm1 = pm1, p
    return pm1

# ── Coulomb amplitude f_C(θ) [fm] ─────────────────────────────────────────────
def coulomb_amp(theta_deg):
    """Standard Rutherford-Coulomb amplitude."""
    th     = np.radians(theta_deg)
    sin_h  = np.sin(th / 2.0)
    sin2   = sin_h * sin_h
    sig0   = sigma_L(0)
    phase  = 2.0 * sig0 - 2.0 * eta * np.log(sin_h)
    return -eta / (2.0 * k * sin2) * np.exp(1j * phase)

# ── Nuclear amplitude f_N(θ) [fm] ─────────────────────────────────────────────
def nuclear_amp(theta_deg, SMatrix):
    """Standard nuclear partial-wave amplitude."""
    cos_th = np.cos(np.radians(theta_deg))
    fN = 0 + 0j
    for L, SL in enumerate(SMatrix):
        sig = sigma_L(L)
        PL  = legendre(L, cos_th)
        fN += (2*L + 1) * np.exp(2j * sig) * (SL - 1.0) * PL
    return fN / (2j * k)

# ── Load S-matrices ────────────────────────────────────────────────────────────
cpp_data  = np.loadtxt('/home/node/working/ptolemy_2019/Cpp_AI/sm148_aa_smat_cpp.txt',  comments='#')
raph_data = np.loadtxt('/home/node/working/ptolemy_2019/Cpp_AI/Raphael_AI/sm148_aa_smat_raphael.txt', comments='#')

S_cpp  = [complex(r[1], r[2]) for r in cpp_data]
S_raph = [complex(r[1], r[2]) for r in raph_data]

print(f"S-matrix lengths: C++={len(S_cpp)}, Raphael={len(S_raph)}")
print(f"  S_cpp[0]  = {S_cpp[0]:.6e}")
print(f"  S_raph[0] = {S_raph[0]:.6e}")
print()

# ── Load C++ computed DCS (for comparison) ────────────────────────────────────
cpp_xsec_data = np.loadtxt('/home/node/working/ptolemy_2019/Cpp_AI/sm148_aa_xsec_cpp.txt', comments='#')
cpp_xsec = {int(round(r[0])): r[1] for r in cpp_xsec_data}

# ── Ptolemy reference ─────────────────────────────────────────────────────────
# From sm148_aa_180_ptolemy.out SIGMA C.M. column [mb/sr]
ptol_ref = {
    10:  1.604e+05,
    20:  2.811e+03,
    30:  1.295e+02,
    40:  2.341e+01,
    60:  4.377e-01,
    90:  3.128e-03,
    120: 1.472e-04,
    150: 2.948e-05,
    180: 4.506e-05,
}

# ── Compute DCS at key angles ─────────────────────────────────────────────────
angles = [5, 10, 20, 30, 40, 60, 90, 120, 150, 180]

# Precompute Coulomb phases for all L (cache)
sigma_cache = [sigma_L(L) for L in range(max(len(S_cpp), len(S_raph)))]
print(f"Coulomb phase σ_0 = {sigma_cache[0]:.6f} rad")
print(f"Coulomb phase σ_1 = {sigma_cache[1]:.6f} rad")
print(f"Coulomb phase σ_5 = {sigma_cache[5]:.6f} rad")
print()

def nuclear_amp_fast(theta_deg, SMatrix):
    cos_th = np.cos(np.radians(theta_deg))
    fN = 0 + 0j
    for L, SL in enumerate(SMatrix):
        sig = sigma_cache[L]
        PL  = legendre(L, cos_th)
        fN += (2*L + 1) * np.exp(2j * sig) * (SL - 1.0) * PL
    return fN / (2j * k)

print(f"{'theta':>6} {'Py_DCS(cpp_S)':>14} {'Py_DCS(raph_S)':>15} "
      f"{'C++_DCS':>12} {'Raphael_DCS':>12} {'Ptolemy':>12} "
      f"{'Py/C++':>8} {'Py/Ptol':>8}")
print("-" * 105)

results = []
for theta in angles:
    fC       = coulomb_amp(theta)
    fN_cpp   = nuclear_amp_fast(theta, S_cpp)
    fN_raph  = nuclear_amp_fast(theta, S_raph)
    
    # DCS in mb/sr: |f|^2 [fm^2] * 10 [fm^2 → mb]
    dcs_py_cpp  = abs(fC + fN_cpp)**2  * 10.0
    dcs_py_raph = abs(fC + fN_raph)**2 * 10.0
    
    dcs_cpp_file  = cpp_xsec.get(theta, float('nan'))
    ptol          = ptol_ref.get(theta, float('nan'))
    
    ratio_py_cpp  = dcs_py_cpp / dcs_cpp_file  if dcs_cpp_file > 0 else float('nan')
    ratio_py_ptol = dcs_py_cpp / ptol          if ptol > 0         else float('nan')
    
    print(f"{theta:6d} {dcs_py_cpp:14.4e} {dcs_py_raph:15.4e} "
          f"{dcs_cpp_file:12.4e} {float('nan'):12.4e} {ptol:12.4e} "
          f"{ratio_py_cpp:8.4f} {ratio_py_ptol:8.4f}")
    
    results.append({
        'theta': theta,
        'dcs_py_cpp': dcs_py_cpp,
        'dcs_py_raph': dcs_py_raph,
        'dcs_cpp': dcs_cpp_file,
        'ptol': ptol,
    })

print()

# ── Also get Raphael DCS file for comparison ───────────────────────────────────
raph_xsec_data = np.loadtxt('/home/node/working/ptolemy_2019/Cpp_AI/Raphael_AI/sm148_aa_xsec_raphael.txt', comments='#')
raph_xsec = {int(round(r[0])): r[1] for r in raph_xsec_data}

print(f"{'theta':>6} {'Py_DCS(cpp_S)':>14} {'Py_DCS(raph_S)':>15} "
      f"{'C++_DCS':>12} {'Raph_DCS':>12} {'Ptolemy':>12} "
      f"{'Py(cpp)/C++':>12} {'Py(raph)/Raph':>14} {'Py/Ptol':>10}")
print("-" * 120)

for r in results:
    theta = r['theta']
    raph_file = raph_xsec.get(theta, float('nan'))
    ratio_py_cpp  = r['dcs_py_cpp']  / r['dcs_cpp']  if r['dcs_cpp']  > 0 else float('nan')
    ratio_py_raph = r['dcs_py_raph'] / raph_file      if raph_file     > 0 else float('nan')
    ratio_py_ptol = r['dcs_py_cpp']  / r['ptol']      if r['ptol']     > 0 else float('nan')
    
    print(f"{theta:6d} {r['dcs_py_cpp']:14.4e} {r['dcs_py_raph']:15.4e} "
          f"{r['dcs_cpp']:12.4e} {raph_file:12.4e} {r['ptol']:12.4e} "
          f"{ratio_py_cpp:12.4f} {ratio_py_raph:14.4f} {ratio_py_ptol:10.4f}")

# ── Diagnose agreement ──────────────────────────────────────────────────────────
print()
print("=" * 60)
print("DIAGNOSIS")
print("=" * 60)

for r in results:
    theta = r['theta']
    raph_file = raph_xsec.get(theta, float('nan'))
    ratio_cpp  = r['dcs_py_cpp']  / r['dcs_cpp']
    ratio_raph = r['dcs_py_raph'] / raph_file if raph_file > 0 else float('nan')
    ratio_ptol = r['dcs_py_cpp']  / r['ptol']
    
    if abs(ratio_cpp - 1) > 0.01:
        print(f"  [BUG] theta={theta}°: Py(cpp_S)/C++DCS = {ratio_cpp:.4f} (>1% off → C++ DCS formula issue)")
    if abs(ratio_raph - 1) > 0.01:
        print(f"  [BUG] theta={theta}°: Py(raph_S)/Raph_DCS = {ratio_raph:.4f} (>1% off → Raphael DCS formula issue)")
    if abs(ratio_ptol - 1) > 0.05:
        print(f"  [DIFF] theta={theta}°: Py(cpp_S)/Ptolemy = {ratio_ptol:.4f} (>5% off vs Ptolemy)")

print()
print(f"S-matrix self-consistency (Py vs C++ DCS using same S-matrix):")
ratios = []
for r in results:
    if r['dcs_cpp'] > 0:
        ratios.append(r['dcs_py_cpp'] / r['dcs_cpp'])
print(f"  mean ratio = {np.mean(ratios):.6f}, std = {np.std(ratios):.6f}")

print(f"\nC++ vs Ptolemy:")
ratios_ptol = []
for r in results:
    if r['ptol'] > 0 and r['dcs_cpp'] > 0:
        ratios_ptol.append(r['dcs_cpp'] / r['ptol'])
print(f"  mean ratio = {np.mean(ratios_ptol):.6f}, std = {np.std(ratios_ptol):.6f}")

print()
print("Done.")
