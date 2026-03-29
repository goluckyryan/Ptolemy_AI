#!/usr/bin/env python3
"""
Plot 2D contour of I_accum magnitude % error between C++ and Fortran
as a function of (Li, Lo).

For each (Li, Lo) pair, we sum |I_accum|^2 over all (JPI, JPO, Lx) channels
in both codes, then compute the % difference of sqrt(sum).
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

def parse_iacc(filename, prefix):
    """Parse I_accum lines from file."""
    data = []
    if prefix == 'CPP':
        pat = re.compile(
            r'CPP_IACC\s+Li=\s*(\d+)\s+JPI=\s*(\d+)\s+Lo=\s*(\d+)\s+JPO=\s*(\d+)\s+Lx=\s*(\d+)\s+'
            r'Ire=\s*([\d.eE+-]+)\s+Iim=\s*([\d.eE+-]+)'
        )
    else:
        pat = re.compile(
            r'FTN_IACC\s+Li=\s*(\d+)\s+JPI=\s*(\d+)\s+Lo=\s*(\d+)\s+JPO=\s*(\d+)\s+Lx=\s*(\d+)\s+'
            r'Ire=\s*([\d.eE+-]+)\s+Iim=\s*([\d.eE+-]+)'
        )
    with open(filename) as f:
        for line in f:
            m = pat.search(line)
            if m:
                Li, JPI, Lo, JPO, Lx = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5))
                Ire, Iim = float(m.group(6)), float(m.group(7))
                data.append((Li, Lo, JPI, JPO, Lx, Ire, Iim))
    return data

# Parse both
cpp_data = parse_iacc('/tmp/cpp_iacc_err.txt', 'CPP')
ftn_data = parse_iacc('/tmp/ftn_iacc_err.txt', 'FTN')

print(f"C++ entries: {len(cpp_data)}")
print(f"Fortran entries: {len(ftn_data)}")

# Build lookup by (Li, Lo, JPI, JPO, Lx) → (Ire, Iim)
cpp_dict = {}
for Li, Lo, JPI, JPO, Lx, Ire, Iim in cpp_data:
    key = (Li, Lo, JPI, JPO, Lx)
    cpp_dict[key] = (Ire, Iim)

ftn_dict = {}
for Li, Lo, JPI, JPO, Lx, Ire, Iim in ftn_data:
    key = (Li, Lo, JPI, JPO, Lx)
    ftn_dict[key] = (Ire, Iim)

# Find all common keys
common_keys = set(cpp_dict.keys()) & set(ftn_dict.keys())
cpp_only = set(cpp_dict.keys()) - set(ftn_dict.keys())
ftn_only = set(ftn_dict.keys()) - set(cpp_dict.keys())
print(f"Common: {len(common_keys)}, C++ only: {len(cpp_only)}, Fortran only: {len(ftn_only)}")

if ftn_only:
    print("Fortran-only entries (first 10):")
    for k in sorted(ftn_only)[:10]:
        print(f"  Li={k[0]} Lo={k[1]} JPI={k[2]} JPO={k[3]} Lx={k[4]}: {ftn_dict[k]}")

# Aggregate by (Li, Lo): sum |I|^2 over (JPI, JPO, Lx)
cpp_LiLo = defaultdict(float)
ftn_LiLo = defaultdict(float)
err_LiLo = defaultdict(list)  # individual magnitude errors

for key in common_keys:
    Li, Lo = key[0], key[1]
    cr, ci = cpp_dict[key]
    fr, fi = ftn_dict[key]
    cpp_mag2 = cr**2 + ci**2
    ftn_mag2 = fr**2 + fi**2
    cpp_LiLo[(Li, Lo)] += cpp_mag2
    ftn_LiLo[(Li, Lo)] += ftn_mag2
    
    # Individual magnitude error
    cpp_mag = np.sqrt(cpp_mag2)
    ftn_mag = np.sqrt(ftn_mag2)
    if ftn_mag > 1e-15:
        err_LiLo[(Li, Lo)].append((cpp_mag - ftn_mag) / ftn_mag * 100)

# Compute % error of RMS magnitude for each (Li, Lo)
all_Li = sorted(set(k[0] for k in cpp_LiLo))
all_Lo = sorted(set(k[1] for k in cpp_LiLo))
print(f"\nLi range: {min(all_Li)}-{max(all_Li)}")
print(f"Lo range: {min(all_Lo)}-{max(all_Lo)}")

# Create 2D grid
Li_max = max(all_Li) + 1
Lo_max = max(all_Lo) + 1
error_grid = np.full((Li_max, Lo_max), np.nan)
abs_error_grid = np.full((Li_max, Lo_max), np.nan)

for (Li, Lo), cpp_sum2 in cpp_LiLo.items():
    ftn_sum2 = ftn_LiLo[(Li, Lo)]
    cpp_rms = np.sqrt(cpp_sum2)
    ftn_rms = np.sqrt(ftn_sum2)
    if ftn_rms > 1e-15:
        pct_err = (cpp_rms - ftn_rms) / ftn_rms * 100
        error_grid[Li, Lo] = pct_err
        abs_error_grid[Li, Lo] = abs(pct_err)

# Print table of errors
print("\n(Li, Lo) → % error (C++ vs Fortran RMS |I_accum|):")
print(f"{'Li':>4} {'Lo':>4} {'%err':>8} {'C++ RMS':>12} {'Ftn RMS':>12}")
for (Li, Lo) in sorted(cpp_LiLo.keys()):
    cpp_rms = np.sqrt(cpp_LiLo[(Li, Lo)])
    ftn_rms = np.sqrt(ftn_LiLo[(Li, Lo)])
    if ftn_rms > 1e-15:
        pct = (cpp_rms - ftn_rms) / ftn_rms * 100
        print(f"{Li:4d} {Lo:4d} {pct:8.2f}% {cpp_rms:12.6f} {ftn_rms:12.6f}")

# Plot 2D contour
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Left: signed error
ax = axes[0]
# Only plot where we have data
valid = ~np.isnan(error_grid)
if valid.any():
    # Use imshow with origin='lower'
    masked = np.ma.masked_invalid(error_grid)
    vmax = max(abs(np.nanmin(error_grid)), abs(np.nanmax(error_grid)))
    vmax = min(vmax, 50)  # cap at 50%
    im = ax.imshow(masked.T, origin='lower', aspect='auto', 
                   cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                   extent=[-0.5, Li_max-0.5, -0.5, Lo_max-0.5])
    plt.colorbar(im, ax=ax, label='% error (C++ − Ftn) / Ftn')
    ax.set_xlabel('Li (incoming)')
    ax.set_ylabel('Lo (outgoing)')
    ax.set_title('I_accum signed % error\n(RMS over JPI, JPO, Lx)')
    ax.set_xlim(0.5, 17.5)
    ax.set_ylim(0.5, 17.5)

# Right: absolute error (log scale)
ax = axes[1]
if valid.any():
    masked_abs = np.ma.masked_invalid(abs_error_grid)
    # Log scale: replace 0 with small number
    log_err = np.where(abs_error_grid > 0, abs_error_grid, 1e-3)
    log_err = np.ma.masked_invalid(log_err)
    im2 = ax.imshow(np.log10(log_err).T, origin='lower', aspect='auto',
                    cmap='hot_r',
                    extent=[-0.5, Li_max-0.5, -0.5, Lo_max-0.5])
    cbar = plt.colorbar(im2, ax=ax, label='log₁₀(|% error|)')
    ax.set_xlabel('Li (incoming)')
    ax.set_ylabel('Lo (outgoing)')
    ax.set_title('I_accum |% error|\n(RMS over JPI, JPO, Lx)')
    ax.set_xlim(0.5, 17.5)
    ax.set_ylim(0.5, 17.5)

plt.suptitle('16O(d,p)17O GS, Elab=20 MeV — I_accum(Li, Lo) C++ vs Fortran', fontsize=14)
plt.tight_layout()
plt.savefig('/home/node/.openclaw/workspace/plots/iacc_2d_error.png', dpi=150)
print("\nPlot saved to plots/iacc_2d_error.png")
