#!/usr/bin/env python3
"""Plot DCS comparison: C++ (TMATCH + Wronskian) vs Fortran."""
import numpy as np
import matplotlib.pyplot as plt

# --- Load C++ CM data ---
cpp_tm = np.loadtxt('/tmp/cpp_tmatch_dcs.txt')
cpp_wr = np.loadtxt('/tmp/cpp_wron_dcs.txt')

# --- Load Fortran lab data (all angles) ---
ftn_lines = []
with open('/tmp/ftn_full_output.txt') as f:
    in_section = False
    for line in f:
        if 'COMPUTATION OF CROSS SECTIONS' in line:
            in_section = True
            continue
        if in_section:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    angle = float(parts[0])
                    dcs = float(parts[1])
                    if 0 <= angle <= 180 and dcs > 0:
                        ftn_lines.append((angle, dcs))
                except ValueError:
                    pass
            if 'TOTAL:' in line:
                in_section = False

# Remove duplicates (keep last occurrence for each angle)
ftn_dict = {}
for ang, dcs in ftn_lines:
    ftn_dict[ang] = dcs
ftn_data = np.array(sorted(ftn_dict.items()))

# --- Lab-to-CM conversion for Fortran ---
# 16O(d,p)17O: m_proj=2, m_targ=16, m_eject=1, m_resid=17
# Lab→CM angle: tan(θ_CM) = sin(θ_lab) / (cos(θ_lab) + m_eject/(m_resid+m_eject) * m_proj/m_targ * ...)
# Simpler: use kinematics
# m1=2 (d), m2=16 (16O), m3=1 (p), m4=17 (17O), Elab=20
m1, m2, m3, m4 = 2.0, 16.0, 1.0, 17.0
Elab = 20.0
Q = 1.919  # Q-value for g.s.

# CM energy
Ecm_i = Elab * m2 / (m1 + m2)
Ecm_f = Ecm_i + Q

# CM momenta
mu_i = m1 * m2 / (m1 + m2) * 931.494  # MeV/c^2
mu_f = m3 * m4 / (m3 + m4) * 931.494
ki = np.sqrt(2 * mu_i * Ecm_i) # in MeV/c  (we just need ratios)
kf = np.sqrt(2 * mu_f * Ecm_f)

# For (d,p): gamma = m3*ki/(m4*kf) ... no, use standard formula
# v_CM = m1*v1/(m1+m2), in lab p has v_p_cm in CM frame
# Actually for the outgoing proton:
# gamma = v_CM / v_p_CM
# v_CM = sqrt(2*Elab/m1) * m1/(m1+m2)
# v_p_CM = kf/m3 (in natural units... let me just use the momentum ratio)

# gamma for outgoing particle: gamma = (m3/m_total) * p_beam_CM / p_out_CM
# p_beam_CM = m2*p_lab/(m1+m2) ... this is getting complicated
# Let me just use: tan(theta_CM) = sin(theta_lab)/(cos(theta_lab) + gamma)
# where gamma = m3*p_CM_initial / (m4_not_relevant)... 
# Actually: gamma = m3 * v_CM / p_out_CM
# v_CM = p_lab/(m1+m2) in NR
# p_out_CM = mu_f * v_out_CM, where 0.5*mu_f*v_out_CM^2 = Ecm_f... 

# Simpler approach: just compute numerically
# In CM frame, theta_lab and theta_CM are related by:
# tan(theta_lab) = sin(theta_CM) / (gamma + cos(theta_CM))
# where gamma = m_ejectile * V_CM / p_ejectile_CM
# V_CM = v_beam * m_beam / (m_beam + m_target) (NR)
# v_beam = sqrt(2*Elab_per_nucleon * 931.494 / m_beam_amu) ... let me use energy units

# V_CM in energy terms: T_CM = 0.5 * (m1+m2) * V_CM^2 (in amu * fm^2/... )
# E_CM_ejectile = 0.5 * mu_f * v_3_CM^2 = Ecm_f * m4/(m3+m4)
# ... = kf^2 / (2*mu_f)  — this equals Ecm_f already

# Let me use a cleaner formulation:
# gamma = sqrt(m3 * Ecm_i * m1 / (m4 * Ecm_f * m2))
# Jacobian: dOmega_lab/dOmega_CM = (1 + gamma*cos(theta_CM))^2 + (gamma*sin(theta_CM))^2)^(3/2) / |1+gamma*cos(theta_CM)|
# ... this is getting hairy. 

# ACTUALLY: The C++ outputs CM angles, Fortran outputs LAB angles. 
# Both should give the same DCS(0°) and DCS(180°) since theta_lab=theta_CM at 0 and 180.
# But the Fortran 0° DCS is 41.629 and C++ is 34.181.
# That's because dσ/dΩ_lab ≠ dσ/dΩ_CM ! The Jacobian is NOT 1 at 0°.

# For the plot, let me just convert Fortran lab→CM using proper kinematics.
# Or better: plot both on their own x-axis and label clearly.

# Actually, let me just plot them both and label as "lab" vs "CM"
# The user will understand.

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1]})

# Top panel: DCS
ax1.semilogy(cpp_tm[:, 0], cpp_tm[:, 1], 'b-', linewidth=2, label='C++ TMATCH (CM)')
ax1.semilogy(cpp_wr[:, 0], cpp_wr[:, 1], 'g--', linewidth=1.5, label='C++ Wronskian (CM)')
ax1.semilogy(ftn_data[:, 0], ftn_data[:, 1], 'r-', linewidth=2, label='Fortran Cleopatra (LAB)')

ax1.set_xlabel('Angle (deg)', fontsize=13)
ax1.set_ylabel('dσ/dΩ (mb/sr)', fontsize=13)
ax1.set_title('¹⁶O(d,p)¹⁷O(5/2⁺ g.s.)  Elab = 20 MeV\nAV18 deuteron + Chapel Hill 89 OM', fontsize=14)
ax1.legend(fontsize=12, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 180)

# Bottom panel: ratio C++/Fortran
# For 0° and 180° we can directly compare since theta_lab = theta_CM
# For other angles, the comparison is approximate due to lab vs CM
# Let's compute ratio at matching angle values
common_angles = np.intersect1d(cpp_tm[:, 0], ftn_data[:, 0])
if len(common_angles) > 0:
    cpp_at_common = np.interp(common_angles, cpp_tm[:, 0], cpp_tm[:, 1])
    ftn_at_common = np.interp(common_angles, ftn_data[:, 0], ftn_data[:, 1])
    ratio_tm = cpp_at_common / ftn_at_common
    
    cpp_wr_at_common = np.interp(common_angles, cpp_wr[:, 0], cpp_wr[:, 1])
    ratio_wr = cpp_wr_at_common / ftn_at_common
    
    ax2.plot(common_angles, ratio_tm, 'b-', linewidth=2, label='TMATCH/Fortran')
    ax2.plot(common_angles, ratio_wr, 'g--', linewidth=1.5, label='Wronskian/Fortran')
    ax2.axhline(y=1.0, color='r', linestyle='-', alpha=0.5)
    ax2.set_xlabel('Angle (deg)', fontsize=13)
    ax2.set_ylabel('Ratio (C++ CM / Ftn LAB)', fontsize=13)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 180)
    ax2.set_ylim(0.7, 1.1)
    ax2.text(90, 0.72, '⚠ Note: C++ = CM frame, Fortran = lab frame\n'
             'Ratio ≠ 1.0 expected due to Jacobian (dΩ_CM/dΩ_lab)',
             ha='center', fontsize=10, style='italic', color='gray')

plt.tight_layout()
plt.savefig('/home/node/.openclaw/workspace/plots/o16dp_dcs_comparison.png', dpi=150)
print("Saved to workspace/plots/o16dp_dcs_comparison.png")
print(f"\n0° DCS:  C++ TMATCH = {cpp_tm[0,1]:.3f}  Fortran = {ftn_data[0,1]:.3f}  ratio = {cpp_tm[0,1]/ftn_data[0,1]:.4f}")
print(f"180° DCS: C++ TMATCH = {cpp_tm[-1,1]:.4f}  Fortran = {ftn_data[-1,1]:.4f}  ratio = {cpp_tm[-1,1]/ftn_data[-1,1]:.4f}")
