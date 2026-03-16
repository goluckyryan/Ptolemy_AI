#!/usr/bin/env python3
"""
148Sm(α,α) at 50 MeV — DCS comparison: C++, Raphael, Ptolemy
Log-scale y-axis, 1°–180°
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re, os

BASE = os.path.dirname(os.path.abspath(__file__))
RAI  = os.path.join(BASE, "Raphael_AI")

# ── C++ DCS ──────────────────────────────────────────────────────────────
cpp = np.loadtxt(os.path.join(BASE, "sm148_aa_xsec_cpp.txt"), comments="#")
cpp_ang, cpp_dcs = cpp[:,0], cpp[:,1]

# ── Raphael DCS ───────────────────────────────────────────────────────────
raph = np.loadtxt(os.path.join(RAI, "sm148_aa_xsec_raphael.txt"), comments="#")
raph_ang, raph_dcs = raph[:,0], raph[:,1]

# ── Ptolemy DCS (from sm148_aa_180_ptolemy.out) ───────────────────────────
ptol_file = os.path.join(RAI, "sm148_aa_180_ptolemy.out")
if not os.path.exists(ptol_file):
    ptol_file = os.path.join(RAI, "sm148_aa_ptolemy.out")

ptol_ang, ptol_dcs = [], []
# format: "   1.00   0.97     1.000        0.1450E+10  ..."
# cols: CM_angle  LAB_angle  sig/Ruth  sig_CM  sig_LAB  Ruth_CM ...
pat = re.compile(r'^\s+(\d+\.\d+)\s+\d+\.\d+\s+\S+\s+([\d.E+\-]+)')
with open(ptol_file) as f:
    for line in f:
        m = pat.match(line)
        if m:
            ang = float(m.group(1))
            dcs = float(m.group(2))
            if dcs > 0:
                ptol_ang.append(ang)
                ptol_dcs.append(dcs)
ptol_ang = np.array(ptol_ang)
ptol_dcs = np.array(ptol_dcs)

print(f"C++:    {len(cpp_ang)} angles, range {cpp_ang[0]:.0f}°–{cpp_ang[-1]:.0f}°")
print(f"Raphael:{len(raph_ang)} angles, range {raph_ang[0]:.0f}°–{raph_ang[-1]:.0f}°")
print(f"Ptolemy:{len(ptol_ang)} angles, range {ptol_ang[0]:.0f}°–{ptol_ang[-1]:.0f}°")

# ── Plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_title(r"$^{148}$Sm($\alpha$,$\alpha$) 50 MeV — Differential Cross Section",
             fontsize=14)

mask_r = raph_dcs > 0
mask_c = cpp_dcs  > 0
mask_p = ptol_dcs > 0

ax.semilogy(raph_ang[mask_r], raph_dcs[mask_r], 'g--', lw=1.8,
            label="Raphael (Python)", zorder=3)
ax.semilogy(ptol_ang[mask_p], ptol_dcs[mask_p], 'ro',  ms=3.5, alpha=0.8,
            label="Ptolemy (Fortran)", zorder=5)
ax.semilogy(cpp_ang[mask_c],  cpp_dcs[mask_c],  'b-',  lw=1.8,
            label="C++ (this work)", zorder=4)

ax.set_xlabel("θ$_{CM}$ (degrees)", fontsize=12)
ax.set_ylabel("dσ/dΩ (mb/sr)", fontsize=12)
ax.set_xlim(0, 180)
ax.legend(fontsize=11)
ax.grid(True, which='both', alpha=0.3)
plt.tight_layout()

out = os.path.join(BASE, "sm148_aa_xsec_all3.png")
plt.savefig(out, dpi=150)
print(f"Saved: {out}")
plt.close()
