#!/usr/bin/env python3
"""
compare_smat_o16dp.py
Compare elastic S-matrix from C++ ElasticSolver vs Ptolemy for both channels
of 16O(d,p)17O at Elab=20 MeV.
"""

import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# -------------------------------------------------------------------------
# Parse Ptolemy .out file for S-matrix
# -------------------------------------------------------------------------
def parse_ptolemy_smat(filename):
    """
    Parse lines like:
      ELASTIC S-MATRIX FOR L =  N, JP =  M/2:  Re +- Im I
    Returns dict keyed by (L, J) -> (Re, Im)
    Ptolemy convention: JP = 2J (integer)
    """
    smat = {}
    # Pattern: L, JP (integer numerator), Re, sign, Im
    pat = re.compile(
        r'ELASTIC S-MATRIX FOR L\s*=\s*(\d+),\s*JP\s*=\s*(-?\d+)/2:\s*'
        r'([+\-]?[\d.Ee+\-]+)\s*\+\s*([+\-]?[\d.Ee+\-]+)\s*I'
    )
    with open(filename) as f:
        for line in f:
            m = pat.search(line)
            if m:
                L    = int(m.group(1))
                twoJ = int(m.group(2))
                Re   = float(m.group(3))
                Im   = float(m.group(4))
                J    = twoJ / 2.0
                smat[(L, J)] = (Re, Im)
    return smat

# -------------------------------------------------------------------------
# Parse C++ .txt file
# -------------------------------------------------------------------------
def parse_cpp_smat(filename):
    """
    Columns: L  J  Re(S)  Im(S)  |S|
    Returns dict keyed by (L, J) -> (Re, Im)
    Skips lines with Re==Im==0 (unphysical entries)
    """
    smat = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            L  = int(parts[0])
            J  = float(parts[1])
            Re = float(parts[2])
            Im = float(parts[3])
            # Skip zero entries (unphysical, e.g. J=L-S for L=0 S=1)
            if Re == 0.0 and Im == 0.0:
                continue
            smat[(L, J)] = (Re, Im)
    return smat

# -------------------------------------------------------------------------
# Compare and print table
# -------------------------------------------------------------------------
def compare_and_print(ptolemy_smat, cpp_smat, channel_name, nprint=15):
    print(f"\n{'='*80}")
    print(f"  {channel_name}")
    print(f"{'='*80}")
    print(f"{'L':>3}  {'J':>5}  {'Pt_Re':>10}  {'Cpp_Re':>10}  {'dRe':>10}  "
          f"{'Pt_Im':>10}  {'Cpp_Im':>10}  {'dIm':>10}  {'|Pt_S|':>8}  {'|Cpp_S|':>8}")
    print("-"*100)

    keys = sorted(ptolemy_smat.keys())
    rows = []
    for (L, J) in keys:
        if (L, J) not in cpp_smat:
            continue
        pt_Re, pt_Im = ptolemy_smat[(L, J)]
        cp_Re, cp_Im = cpp_smat[(L, J)]
        dRe = cp_Re - pt_Re
        dIm = cp_Im - pt_Im
        pt_abs = np.sqrt(pt_Re**2 + pt_Im**2)
        cp_abs = np.sqrt(cp_Re**2 + cp_Im**2)
        rows.append((L, J, pt_Re, cp_Re, dRe, pt_Im, cp_Im, dIm, pt_abs, cp_abs))

    # Print first nprint rows
    for row in rows[:nprint]:
        L, J, pt_Re, cp_Re, dRe, pt_Im, cp_Im, dIm, pt_abs, cp_abs = row
        print(f"{L:>3}  {J:>5.1f}  {pt_Re:>10.5f}  {cp_Re:>10.5f}  {dRe:>+10.5f}  "
              f"{pt_Im:>10.5f}  {cp_Im:>10.5f}  {dIm:>+10.5f}  "
              f"{pt_abs:>8.5f}  {cp_abs:>8.5f}")

    if rows:
        dRes = np.array([r[4] for r in rows])
        dIms = np.array([r[7] for r in rows])
        print(f"\n  Max |ΔRe(S)| = {np.max(np.abs(dRes)):.6f}")
        print(f"  Max |ΔIm(S)| = {np.max(np.abs(dIms)):.6f}")
        print(f"  RMS  ΔRe(S)  = {np.sqrt(np.mean(dRes**2)):.6f}")
        print(f"  RMS  ΔIm(S)  = {np.sqrt(np.mean(dIms**2)):.6f}")
        
        # Check <1% criterion (against |S| ~ 1 scale)
        max_diff = max(np.max(np.abs(dRes)), np.max(np.abs(dIms)))
        if max_diff < 0.01:
            print(f"  ✓ Agreement < 1% (max diff = {max_diff:.4f})")
        else:
            print(f"  ✗ Max difference = {max_diff:.4f} > 1%")
    return rows

# -------------------------------------------------------------------------
# Generate plots
# -------------------------------------------------------------------------
def make_plots(ptolemy_smat, cpp_smat, channel, png_dir):
    os.makedirs(png_dir, exist_ok=True)
    keys = sorted(set(ptolemy_smat.keys()) & set(cpp_smat.keys()))

    # Separate J+/J- series for spin-1/2, or just all for spin-1
    # Determine spin from J values
    Js = sorted(set(J for (L, J) in keys))
    
    # Build arrays: one per J = L ± S
    # For each (L, J) pair, collect data
    data_pt = {k: ptolemy_smat[k] for k in keys}
    data_cp = {k: cpp_smat[k]     for k in keys}

    Ls  = [k[0] for k in keys]
    Js_ = [k[1] for k in keys]
    pt_Re = [data_pt[k][0] for k in keys]
    pt_Im = [data_pt[k][1] for k in keys]
    cp_Re = [data_cp[k][0] for k in keys]
    cp_Im = [data_cp[k][1] for k in keys]
    dRe   = [cp_Re[i]-pt_Re[i] for i in range(len(keys))]
    dIm   = [cp_Im[i]-pt_Im[i] for i in range(len(keys))]

    # Use unique L values for the x-axis label
    uniq_L = sorted(set(Ls))

    # For clarity, create an index for each (L, J) pair
    idx = list(range(len(keys)))
    xlabels = [f"L={k[0]}\nJ={k[1]:.1f}" for k in keys]

    # Plot 1: Re and Im S-matrix comparison
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    axes[0].plot(idx, pt_Re, 'b-o', ms=4, label='Ptolemy Re(S)')
    axes[0].plot(idx, cp_Re, 'r--s', ms=4, label='C++ Re(S)')
    axes[0].set_ylabel('Re(S)')
    axes[0].legend(fontsize=9)
    axes[0].set_title(f'{channel}: S-matrix comparison (Re)')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(idx, pt_Im, 'b-o', ms=4, label='Ptolemy Im(S)')
    axes[1].plot(idx, cp_Im, 'r--s', ms=4, label='C++ Im(S)')
    axes[1].set_ylabel('Im(S)')
    axes[1].legend(fontsize=9)
    axes[1].set_title(f'{channel}: S-matrix comparison (Im)')
    axes[1].grid(True, alpha=0.3)

    plt.xticks(idx[::3], xlabels[::3], fontsize=6, rotation=45)
    plt.tight_layout()
    fname = os.path.join(png_dir, f'{channel}_smat_compare.png')
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")

    # Plot 2: differences
    fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
    axes[0].semilogy(idx, np.abs(dRe), 'b-o', ms=4, label='|ΔRe(S)|')
    axes[0].axhline(0.01, color='gray', ls='--', label='1% level')
    axes[0].set_ylabel('|ΔRe(S)|')
    axes[0].legend(fontsize=9)
    axes[0].set_title(f'{channel}: S-matrix differences')
    axes[0].grid(True, alpha=0.3)

    axes[1].semilogy(idx, np.abs(dIm), 'r-s', ms=4, label='|ΔIm(S)|')
    axes[1].axhline(0.01, color='gray', ls='--', label='1% level')
    axes[1].set_ylabel('|ΔIm(S)|')
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)

    plt.xticks(idx[::3], xlabels[::3], fontsize=6, rotation=45)
    plt.tight_layout()
    fname2 = os.path.join(png_dir, f'{channel}_smat_diff.png')
    plt.savefig(fname2, dpi=150)
    plt.close()
    print(f"  Saved {fname2}")

# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
if __name__ == '__main__':
    base = 'Raphael_AI'
    png  = 'png'

    print("Loading Ptolemy outputs...")
    d16O_pt = parse_ptolemy_smat(f'{base}/d16O_elastic.out')
    p17O_pt = parse_ptolemy_smat(f'{base}/p17O_elastic.out')
    print(f"  d+16O: {len(d16O_pt)} entries")
    print(f"  p+17O: {len(p17O_pt)} entries")

    print("Loading C++ outputs...")
    d16O_cp = parse_cpp_smat(f'{base}/d16O_smat_cpp.txt')
    p17O_cp = parse_cpp_smat(f'{base}/p17O_smat_cpp.txt')
    print(f"  d+16O: {len(d16O_cp)} entries")
    print(f"  p+17O: {len(p17O_cp)} entries")

    # Print comparison tables
    rows_d = compare_and_print(d16O_pt, d16O_cp, "d + 16O (Elab=20 MeV, An-Cai 2006)")
    rows_p = compare_and_print(p17O_pt, p17O_cp, "p + 17O (Elab=23.21 MeV, Koning-Delaroche)")

    # Generate plots
    print("\nGenerating plots...")
    make_plots(d16O_pt, d16O_cp, 'd16O', png)
    make_plots(p17O_pt, p17O_cp, 'p17O', png)

    print("\nDone.")
