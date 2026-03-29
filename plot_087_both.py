#!/usr/bin/env python3
"""Plot 0.871 MeV state DCS at 10 and 20 MeV: C++ vs Fortran, CM frame."""
import numpy as np
import matplotlib.pyplot as plt

# 20 MeV
cpp_tm_20 = np.loadtxt('/tmp/cpp_087_tmatch.txt')
cpp_wr_20 = np.loadtxt('/tmp/cpp_087_wron.txt')
ftn_20 = np.loadtxt('/tmp/ftn_087_cm.txt')

# 10 MeV
cpp_tm_10 = np.loadtxt('/tmp/cpp_087_10_tmatch.txt')
cpp_wr_10 = np.loadtxt('/tmp/cpp_087_10_wron.txt')
ftn_10 = np.loadtxt('/tmp/ftn_087_10_cm.txt')

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for col, (E, cpp_tm, cpp_wr, ftn) in enumerate([
    (10, cpp_tm_10, cpp_wr_10, ftn_10),
    (20, cpp_tm_20, cpp_wr_20, ftn_20),
]):
    err_tm = (cpp_tm[:, 1] / ftn[:, 1] - 1) * 100
    err_wr = (cpp_wr[:, 1] / ftn[:, 1] - 1) * 100
    
    print(f"\nElab = {E} MeV:")
    print(f"  TMATCH:    0° = {cpp_tm[0,1]:.3f} vs {ftn[0,1]:.3f}  ({err_tm[0]:.2f}%)")
    print(f"  Wronskian: 0° = {cpp_wr[0,1]:.3f} vs {ftn[0,1]:.3f}  ({err_wr[0]:.2f}%)")
    print(f"  TMATCH:    mean |err| = {np.mean(np.abs(err_tm)):.2f}%  max = {np.max(np.abs(err_tm)):.2f}% at {ftn[np.argmax(np.abs(err_tm)),0]:.0f}°")
    print(f"  Wronskian: mean |err| = {np.mean(np.abs(err_wr)):.2f}%  max = {np.max(np.abs(err_wr)):.2f}% at {ftn[np.argmax(np.abs(err_wr)),0]:.0f}°")
    
    ax = axes[0, col]
    ax.semilogy(ftn[:, 0], ftn[:, 1], 'r-', lw=2, label='Fortran')
    ax.semilogy(cpp_tm[:, 0], cpp_tm[:, 1], 'b--', lw=2, label='C++ TMATCH')
    ax.semilogy(cpp_wr[:, 0], cpp_wr[:, 1], 'g:', lw=2, label='C++ Wronskian')
    ax.set_title(f'Elab = {E} MeV', fontsize=13)
    ax.set_ylabel('dσ/dΩ_CM (mb/sr)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 180)
    
    ax2 = axes[1, col]
    ax2.plot(ftn[:, 0], err_tm, 'b-', lw=1.5, 
             label=f'TMATCH ({np.mean(np.abs(err_tm)):.1f}%)')
    ax2.plot(ftn[:, 0], err_wr, 'g-', lw=1.5,
             label=f'Wronskian ({np.mean(np.abs(err_wr)):.1f}%)')
    ax2.axhline(0, color='r', alpha=0.5)
    ax2.fill_between(ftn[:, 0], -5, 5, alpha=0.08, color='green')
    ax2.set_ylabel('% Error', fontsize=12)
    ax2.set_xlabel('θ_CM (deg)', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 180)
    ax2.set_ylim(-25, 25)

fig.suptitle('¹⁶O(d,p)¹⁷O(1/2⁺ 0.871 MeV) — CM frame\nAV18 deuteron + CH89 OM', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('/home/node/.openclaw/workspace/plots/o16dp_087_10v20.png', dpi=150, bbox_inches='tight')
print("\nSaved.")
