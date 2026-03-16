#!/usr/bin/env python3
"""
Compare 148Sm(a,a) at 50 MeV: Raphael Python vs Ptolemy Fortran
S-matrix (Re, Im, |S|) and dσ/dΩ from 1° to 180° (1° step).

Ptolemy reference: sm148_aa_ptolemy.out
Potential: V=65.5, W=29.8, R0=RI0=1.427, A=AI=0.671 (volume WS), RC0=1.4
"""

import sys
import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- add Raphael to path ----
RAPHAEL = "/home/node/working/PtolemyGUI/Raphael"
sys.path.insert(0, RAPHAEL)

from distortedWave import DistortedWave
from solveSE import CoulombPotential, WoodsSaxonPot

# ============================================================
# 1.  Set up Raphael for 148Sm(a,a) at 50 MeV
# ============================================================
# Use 'a' shorthand (same as 4He in IAEA database)
# Potential: combined complex WS + Coulomb (Raphael convention: useBothMass=False)
dw = DistortedWave("148Sm", "a", 50)
dw.ClearPotential()
# Ptolemy convention for alpha: R = r0 * (At^{1/3} + Ap^{1/3}) — useBothMass=True
# Verified: R0 = 1.427*(148^{1/3}+4^{1/3}) = 9.8134 fm  (matches Ptolemy output)
#            Rc = 1.400*(148^{1/3}+4^{1/3}) = 9.6278 fm  (matches Ptolemy output)
dw.AddPotential(WoodsSaxonPot(-65.5 - 29.8j, 1.427, 0.671), True)
dw.AddPotential(CoulombPotential(1.4), True)

print("=== Kinematics ===")
dw.PrintInput()
print()

# ============================================================
# 2.  Solve S-matrix (alpha is spin-0, so S=0, only J=L)
# ============================================================
# Grid: dr=0.05 fm, rMax=20 fm — matches Ptolemy's ASYMPTOPIA=20 fm exactly.
# Matching at 20 fm is correct: nuclear potential is negligible there.
# Larger rMax causes RK4 phase accumulation errors for high-L near-unity S values.
dw.SetRange(0.0, 0.05, 400)
# Optimal grid: dr=0.05 fm, rMax=25 fm — matches nuclear potential range;
# matching at larger r contaminates S-matrix via mpmath Coulomb function errors at large rho.
dw.SetRange(0.0, 0.05, 500)   # 500 * 0.05 = 25 fm
print("Calculating S-matrix...")
dw.CalScatteringMatrix(maxL=90)

# ============================================================
# 3.  Parse Ptolemy S-matrix output
# ============================================================
PTOL_OUT = os.path.join(os.path.dirname(__file__), "sm148_aa_180_ptolemy.out")

ptol_smat = {}   # L -> (Re, Im)
ptol_xsec = {}   # theta_cm -> sigma_cm (mb/sr)

with open(PTOL_OUT) as f:
    lines = f.readlines()

in_smat = False
in_xsec = False
for line in lines:
    if "L    L'   LX" in line:
        in_smat = True
        in_xsec = False
        continue
    if "ANGLE" in line and "SIGMA" in line:
        in_smat = False
        in_xsec = True
        continue
    if "TOTAL REACTION" in line or "INPUT... END" in line:
        in_smat = False
        in_xsec = False

    if in_smat:
        # "  LL   LL    0     Re +- Im I   |S|  phase  frac"
        m = re.match(r'\s+(\d+)\s+\d+\s+0\s+([-+]?\d+\.\d+E[+-]\d+)\s+([\+-]+)([\d.E+-]+)\s+I', line)
        if m:
            L   = int(m.group(1))
            re_ = float(m.group(2))
            sign = -1.0 if '+-' in m.group(3) else 1.0
            im_ = sign * float(m.group(4))
            ptol_smat[L] = (re_, im_)

    if in_xsec:
        m = re.match(r'\s+([\d.]+)\s+[\d.]+\s+([\d.E+-]+)\s+([\d.E+-]+)', line)
        if m:
            theta = float(m.group(1))
            sigma_cm = float(m.group(3))   # absolute dσ/dΩ in mb/sr (C.M.)
            ptol_xsec[theta] = sigma_cm

# ============================================================
# 4.  Print S-matrix comparison table
# ============================================================
print("\n=== S-matrix comparison (Ptolemy vs Raphael) ===")
print(f"{'L':>4}  {'Ptol Re':>12} {'Ptol Im':>12}  {'Raph Re':>12} {'Raph Im':>12}  {'|S|Ptol':>9} {'|S|Raph':>9}")
print("-" * 85)

for L in sorted(ptol_smat.keys()):
    pr, pi = ptol_smat[L]
    if L < len(dw.ScatMatrix):
        S_r = dw.ScatMatrix[L][0]
        rr, ri = np.real(S_r), np.imag(S_r)
        mag_p = np.sqrt(pr**2 + pi**2)
        mag_r = abs(S_r)
        print(f"{L:4d}  {pr:12.6f} {pi:12.6f}  {rr:12.6f} {ri:12.6f}  {mag_p:9.5f} {mag_r:9.5f}")
    else:
        print(f"{L:4d}  {pr:12.6f} {pi:12.6f}  {'---':>12} {'---':>12}")

# ============================================================
# 5.  Compute Raphael DCS 1°–180° and compare vs Ptolemy 1°–100°
# ============================================================
print("\n=== Angular distribution comparison (0°–180°) ===")
print(f"{'θ_cm':>7}  {'Ptol σ (mb/sr)':>16}  {'Raph σ (mb/sr)':>16}  {'Ratio R/P':>10}")
print("-" * 58)

def dcs_nuclear_only(dw, theta_deg):
    """Compute elastic DCS: f = f_C (analytic) + f_N (plain partial wave sum).
    f_N = (1/2ik) Σ_L (2L+1) e^{2iσ_L} (S_L - 1) P_L(cosθ)
    (S_L - 1) → 0 rapidly for L >> Lcrit so the sum is finite.
    Returns DCS in fm^2."""
    costh = np.cos(np.radians(theta_deg))
    k = dw.k
    fc = dw.CoulombScatterintAmp(theta_deg)
    fn = 0.0 + 0j
    maxL = len(dw.ScatMatrix) - 1
    for L in range(0, maxL + 1):
        SL = dw.ScatMatrix[L][0]
        sigma_L = dw.CoulombPhaseShift(L)
        PL = float(np.polynomial.legendre.legval(costh, [0]*L + [1]))
        fn += (2*L + 1) * np.exp(2j * sigma_L) * (SL - 1) * PL
    fn /= (2j * k)
    return abs(fc + fn)**2

theta_list = list(range(0, 181))
raph_xsec = {}

for theta in theta_list:
    if theta == 0:
        raph_xsec[0] = float('nan')  # Rutherford diverges at 0°
        continue
    raph_xsec[theta] = dcs_nuclear_only(dw, theta) * 10.0  # fm^2 -> mb/sr

for theta in range(1, 181):
    ps = ptol_xsec.get(float(theta), None)
    rs = raph_xsec.get(theta, None)
    if ps is None or rs is None:
        continue
    ratio = rs / ps if ps != 0 else float('nan')
    flag = "  ← !" if abs(ratio - 1) > 0.05 else ""
    print(f"{theta:7.1f}  {ps:16.4e}  {rs:16.4e}  {ratio:10.5f}{flag}")

# ============================================================
# 6.  Plots
# ============================================================
# 6a. S-matrix Re/Im vs L
L_arr = list(range(len(dw.ScatMatrix)))
S_re  = [np.real(dw.ScatMatrix[L][0]) for L in L_arr]
S_im  = [np.imag(dw.ScatMatrix[L][0]) for L in L_arr]
ptol_L  = sorted(ptol_smat.keys())
ptol_re = [ptol_smat[L][0] for L in ptol_L]
ptol_im = [ptol_smat[L][1] for L in ptol_L]

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
axes[0].plot(L_arr,  S_re,    'b-o',  ms=4, lw=1, label='Raphael Re(S)')
axes[0].plot(ptol_L, ptol_re, 'r--s', ms=5, lw=1, label='Ptolemy Re(S)')
axes[0].set_xlabel('L'); axes[0].set_ylabel('Re(S)')
axes[0].set_title('148Sm(α,α) 50 MeV — S-matrix Real Part')
axes[0].legend(); axes[0].grid(True)

axes[1].plot(L_arr,  S_im,    'b-o',  ms=4, lw=1, label='Raphael Im(S)')
axes[1].plot(ptol_L, ptol_im, 'r--s', ms=5, lw=1, label='Ptolemy Im(S)')
axes[1].set_xlabel('L'); axes[1].set_ylabel('Im(S)')
axes[1].set_title('148Sm(α,α) 50 MeV — S-matrix Imaginary Part')
axes[1].legend(); axes[1].grid(True)

plt.tight_layout()
smat_fig = os.path.join(os.path.dirname(__file__), "sm148_aa_smat.png")
plt.savefig(smat_fig, dpi=120)
plt.close()
print(f"\nSaved S-matrix plot: {smat_fig}")

# 6b. DCS — absolute and σ/Rutherford
ruth = {}
for theta in theta_list:
    ruth[theta] = float(np.real(dw.RutherFord(theta))) * 10.0  # fm^2 → mb/sr

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Absolute DCS (1–180°, both Ptolemy & Raphael)
theta_plot = list(range(1, 181))
ptol_vals = [ptol_xsec.get(float(t), np.nan) for t in theta_plot]
raph_vals_plot = [raph_xsec.get(t, np.nan) for t in theta_plot]

axes[0].semilogy(theta_plot, ptol_vals,      'r--', ms=3, lw=1, label='Ptolemy')
axes[0].semilogy(theta_plot, raph_vals_plot, 'b-',  ms=3, lw=1, label='Raphael')
axes[0].set_xlabel('θ_cm (deg)'); axes[0].set_ylabel('dσ/dΩ (mb/sr)')
axes[0].set_title('148Sm(α,α) 50 MeV — Absolute DCS (1°–180°)')
axes[0].legend(); axes[0].grid(True, which='both', alpha=0.4)

# σ/Rutherford (both 1–180°)
ruth = {t: float(np.real(dw.RutherFord(t))) * 10.0 for t in theta_plot}
raph_ratio_full = [raph_xsec.get(t, np.nan) / ruth.get(t, np.nan) for t in theta_plot]
ptol_ratio      = [ptol_xsec.get(float(t), np.nan) / ruth.get(t, np.nan) for t in theta_plot]

axes[1].semilogy(theta_plot, raph_ratio_full, 'b-',  lw=1.5, label='Raphael σ/Ruth')
axes[1].semilogy(theta_plot, ptol_ratio,      'r--', lw=1.5, label='Ptolemy σ/Ruth')
axes[1].set_xlabel('θ_cm (deg)'); axes[1].set_ylabel('σ/σ_Rutherford')
axes[1].set_title('148Sm(α,α) 50 MeV — σ/Rutherford (0°–180°)')
axes[1].legend(); axes[1].grid(True, which='both', alpha=0.4)

plt.tight_layout()
xsec_fig = os.path.join(os.path.dirname(__file__), "sm148_aa_xsec.png")
plt.savefig(xsec_fig, dpi=120)
plt.close()
print(f"Saved DCS plot: {xsec_fig}")
print("\nDone.")
