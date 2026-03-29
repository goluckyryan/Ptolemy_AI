#!/usr/bin/env python3
"""Plot Ryan's 0.87 MeV state: C++ vs Fortran, CM frame."""
import numpy as np
import matplotlib.pyplot as plt

cpp = np.loadtxt('/tmp/cpp_ryan_087.txt')
ftn = np.loadtxt('/tmp/ftn_ryan_087_cm.txt')

err = (cpp[:, 1] / ftn[:, 1] - 1) * 100

print(f"0° DCS: C++={cpp[0,1]:.3f}  Ftn={ftn[0,1]:.3f}  err={err[0]:.2f}%")
print(f"Mean |err| = {np.mean(np.abs(err)):.2f}%  Max = {np.max(np.abs(err)):.2f}% at {ftn[np.argmax(np.abs(err)),0]:.1f}°")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1]})

ax1.semilogy(ftn[:, 0], ftn[:, 1], 'r-', lw=2, label='Fortran Cleopatra')
ax1.semilogy(cpp[:, 0], cpp[:, 1], 'b--', lw=2, label='C++ (Wronskian)')
ax1.set_xlabel('θ_CM (deg)', fontsize=13)
ax1.set_ylabel('dσ/dΩ_CM (mb/sr)', fontsize=13)
ax1.set_title("¹⁶O(d,p)¹⁷O(1/2⁺ 0.87 MeV)  Elab = 20 MeV — CM frame\n"
              "Ryan's input: r0_T=1.25, VSO_T=6, asymptopia=50, lmax=30\n"
              "An&Cai (d) + KD03 (p)", fontsize=13)
ax1.legend(fontsize=12)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 180)

ax2.plot(ftn[:, 0], err, 'b-', lw=1.5, 
         label=f'mean |err|={np.mean(np.abs(err)):.1f}%')
ax2.axhline(0, color='r', alpha=0.5)
ax2.fill_between(ftn[:, 0], -5, 5, alpha=0.08, color='green')
ax2.set_xlabel('θ_CM (deg)', fontsize=13)
ax2.set_ylabel('% Error vs Fortran', fontsize=13)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 180)
ax2.set_ylim(-40, 40)

plt.tight_layout()
plt.savefig('/home/node/.openclaw/workspace/plots/o16dp_ryan_087_comparison.png', dpi=150)
print("Saved.")
