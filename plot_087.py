#!/usr/bin/env python3
"""Plot 0.871 MeV state DCS: C++ vs Fortran, CM frame, 0-180 deg."""
import numpy as np
import matplotlib.pyplot as plt

cpp_tm = np.loadtxt('/tmp/cpp_087_tmatch.txt')
cpp_wr = np.loadtxt('/tmp/cpp_087_wron.txt')
ftn = np.loadtxt('/tmp/ftn_087_cm.txt')

# % error
err_tm = (cpp_tm[:, 1] / ftn[:, 1] - 1) * 100
err_wr = (cpp_wr[:, 1] / ftn[:, 1] - 1) * 100

print(f"TMATCH:    0° DCS = {cpp_tm[0,1]:.3f}  Ftn = {ftn[0,1]:.3f}  err = {err_tm[0]:.2f}%")
print(f"Wronskian: 0° DCS = {cpp_wr[0,1]:.3f}  Ftn = {ftn[0,1]:.3f}  err = {err_wr[0]:.2f}%")
print(f"\nTMATCH:    mean |err| = {np.mean(np.abs(err_tm)):.2f}%  max = {np.max(np.abs(err_tm)):.2f}% at {cpp_tm[np.argmax(np.abs(err_tm)),0]:.0f}°")
print(f"Wronskian: mean |err| = {np.mean(np.abs(err_wr)):.2f}%  max = {np.max(np.abs(err_wr)):.2f}% at {cpp_wr[np.argmax(np.abs(err_wr)),0]:.0f}°")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1]})

ax1.semilogy(ftn[:, 0], ftn[:, 1], 'r-', linewidth=2, label='Fortran Cleopatra')
ax1.semilogy(cpp_tm[:, 0], cpp_tm[:, 1], 'b--', linewidth=2, label='C++ TMATCH')
ax1.semilogy(cpp_wr[:, 0], cpp_wr[:, 1], 'g:', linewidth=2, label='C++ Wronskian')

ax1.set_xlabel('θ_CM (deg)', fontsize=13)
ax1.set_ylabel('dσ/dΩ_CM (mb/sr)', fontsize=13)
ax1.set_title('¹⁶O(d,p)¹⁷O(1/2⁺ 0.871 MeV)  Elab = 20 MeV — CM frame\nAV18 deuteron | An&Cai (d) + KD03 (p)', fontsize=14)
ax1.legend(fontsize=12, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 180)

ax2.plot(ftn[:, 0], err_tm, 'b-', linewidth=1.5, label=f'TMATCH (mean |err|={np.mean(np.abs(err_tm)):.1f}%)')
ax2.plot(ftn[:, 0], err_wr, 'g-', linewidth=1.5, label=f'Wronskian (mean |err|={np.mean(np.abs(err_wr)):.1f}%)')
ax2.axhline(0, color='r', alpha=0.5)
ax2.fill_between(ftn[:, 0], -5, 5, alpha=0.08, color='green')
ax2.set_xlabel('θ_CM (deg)', fontsize=13)
ax2.set_ylabel('% Error vs Fortran', fontsize=13)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 180)
ax2.set_ylim(-25, 25)

plt.tight_layout()
plt.savefig('/home/node/.openclaw/workspace/plots/o16dp_087_dcs_cm.png', dpi=150)
print("\nSaved plot.")
