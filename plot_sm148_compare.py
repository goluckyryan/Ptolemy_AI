#!/usr/bin/env python3
"""
Comparison plots: 148Sm(α,α) at 50 MeV
S-matrix (Re/Im vs L) and DCS (mb/sr vs angle)
Three datasets: C++, Raphael, Ptolemy
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import os

BASEDIR = os.path.dirname(os.path.abspath(__file__))
RAI     = os.path.join(BASEDIR, "Raphael_AI")

# ── Load Raphael S-matrix ──────────────────────────────────────────────────
raph_smat = np.loadtxt(os.path.join(RAI, "sm148_aa_smat_raphael.txt"), comments="#")
raph_L   = raph_smat[:, 0]
raph_re  = raph_smat[:, 1]
raph_im  = raph_smat[:, 2]

# ── Load Raphael DCS ───────────────────────────────────────────────────────
raph_xsec = np.loadtxt(os.path.join(RAI, "sm148_aa_xsec_raphael.txt"), comments="#")
raph_ang  = raph_xsec[:, 0]
raph_dcs  = raph_xsec[:, 1]

# ── Load C++ S-matrix ──────────────────────────────────────────────────────
cpp_smat = np.loadtxt(os.path.join(BASEDIR, "sm148_aa_smat_cpp.txt"), comments="#")
cpp_L   = cpp_smat[:, 0]
cpp_re  = cpp_smat[:, 1]
cpp_im  = cpp_smat[:, 2]

# ── Load C++ DCS ───────────────────────────────────────────────────────────
cpp_xsec = np.loadtxt(os.path.join(BASEDIR, "sm148_aa_xsec_cpp.txt"), comments="#")
cpp_ang  = cpp_xsec[:, 0]
cpp_dcs  = cpp_xsec[:, 1]

# ── Parse Ptolemy output for S-matrix ────────────────────────────────────
ptol_L, ptol_re, ptol_im = [], [], []
pto_file = os.path.join(RAI, "sm148_aa_180_ptolemy.out")
if not os.path.exists(pto_file):
    pto_file = os.path.join(RAI, "sm148_aa_ptolemy.out")

smat_pattern = re.compile(
    r'^\s*(\d+)\s+\d+\s+\d+\s+([-+]?\d+\.\d+E[+-]\d+)\s+\+?(-?\d+\.\d+E[+-]\d+)\s+I'
)
with open(pto_file) as f:
    for line in f:
        line = line.rstrip()
        m = smat_pattern.match(line)
        if m:
            L   = int(m.group(1))
            re_ = float(m.group(2))
            im_ = float(m.group(3))
            ptol_L.append(L)
            ptol_re.append(re_)
            ptol_im.append(im_)

ptol_L  = np.array(ptol_L)
ptol_re = np.array(ptol_re)
ptol_im = np.array(ptol_im)

# ── Parse Ptolemy DCS ─────────────────────────────────────────────────────
pto_xsec_file = os.path.join(RAI, "sm148_aa_180_ptolemy.out")
if not os.path.exists(pto_xsec_file):
    pto_xsec_file = os.path.join(RAI, "sm148_aa_ptolemy.out")

ptol_ang_list, ptol_dcs_list = [], []
# Ptolemy DCS format: "   1.00   0.97     1.000        0.1450E+10 ..."
# columns: CM_angle  LAB_angle  sigma/Ruth  sigma_CM  sigma_LAB  Ruth_CM  ...
# We want CM angle (col 0) and sigma_CM (col 3)
dcs_pattern = re.compile(r'^\s+(\d+\.\d+)\s+\d+\.\d+\s+[-\d.E+]+\s+([-\d.E+]+)')
with open(pto_xsec_file) as f:
    for line in f:
        m = dcs_pattern.match(line)
        if m:
            ang = float(m.group(1))
            dcs = float(m.group(2))
            if ang > 0 and dcs > 0:
                ptol_ang_list.append(ang)
                ptol_dcs_list.append(dcs)

ptol_ang = np.array(ptol_ang_list)
ptol_dcs = np.array(ptol_dcs_list)

print(f"Raphael smat: {len(raph_L)} L values, dcs: {len(raph_ang)} angles")
print(f"C++     smat: {len(cpp_L)} L values, dcs: {len(cpp_ang)} angles")
print(f"Ptolemy smat: {len(ptol_L)} L values, dcs: {len(ptol_ang)} angles")

# ══════════════════════════════════════════════════════════════════════════
# PLOT 1 — S-matrix comparison
# ══════════════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
fig.suptitle(r"$^{148}$Sm($\alpha$,$\alpha$) at 50 MeV — S-matrix", fontsize=13)

# Re(S)
ax1.plot(raph_L,  raph_re,  'g--',  lw=1.5, label="Raphael", zorder=3)
ax1.plot(cpp_L,   cpp_re,   'b-',   lw=1.8, label="C++",     zorder=4)
if len(ptol_L):
    ax1.plot(ptol_L,  ptol_re,  'ro',   ms=4,   label="Ptolemy", zorder=5, alpha=0.8)
ax1.axhline(0, color='k', lw=0.5, ls=':')
ax1.axhline(1, color='k', lw=0.5, ls=':')
ax1.set_ylabel("Re(S)", fontsize=11)
ax1.legend(fontsize=10)
ax1.set_ylim(-0.15, 1.1)
ax1.grid(True, alpha=0.3)

# Im(S)
ax2.plot(raph_L,  raph_im,  'g--',  lw=1.5, label="Raphael", zorder=3)
ax2.plot(cpp_L,   cpp_im,   'b-',   lw=1.8, label="C++",     zorder=4)
if len(ptol_L):
    ax2.plot(ptol_L,  ptol_im,  'ro',   ms=4,   label="Ptolemy", zorder=5, alpha=0.8)
ax2.axhline(0, color='k', lw=0.5, ls=':')
ax2.set_ylabel("Im(S)", fontsize=11)
ax2.set_xlabel("L (angular momentum)", fontsize=11)
ax2.legend(fontsize=10)
ax2.set_ylim(-0.35, 0.35)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
out1 = os.path.join(BASEDIR, "sm148_aa_smat_compare.png")
plt.savefig(out1, dpi=150)
print(f"Saved: {out1}")
plt.close()

# ══════════════════════════════════════════════════════════════════════════
# PLOT 2 — DCS comparison (log scale)
# ══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(9, 6))
ax.set_title(r"$^{148}$Sm($\alpha$,$\alpha$) at 50 MeV — Differential Cross Section", fontsize=13)

# mask zeros/negatives for log
mask_r = raph_dcs > 0
mask_c = cpp_dcs  > 0

ax.semilogy(raph_ang[mask_r], raph_dcs[mask_r], 'g--', lw=1.5, label="Raphael", zorder=3)
ax.semilogy(cpp_ang[mask_c],  cpp_dcs[mask_c],  'b-',  lw=1.8, label="C++",     zorder=4)
if len(ptol_ang):
    mask_p = ptol_dcs > 0
    ax.semilogy(ptol_ang[mask_p], ptol_dcs[mask_p], 'ro', ms=4, label="Ptolemy", zorder=5, alpha=0.8)

ax.set_xlabel("θ (degrees)", fontsize=11)
ax.set_ylabel("dσ/dΩ (mb/sr)", fontsize=11)
ax.set_xlim(0, 180)
ax.legend(fontsize=10)
ax.grid(True, which='both', alpha=0.3)

plt.tight_layout()
out2 = os.path.join(BASEDIR, "sm148_aa_xsec_compare.png")
plt.savefig(out2, dpi=150)
print(f"Saved: {out2}")
plt.close()
