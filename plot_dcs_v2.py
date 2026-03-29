#!/usr/bin/env python3
"""Plot DCS: C++ CM→lab converted vs Fortran lab. Proper Jacobian."""
import numpy as np
import matplotlib.pyplot as plt

# Masses in amu
m1 = 2.01410178   # deuteron
m2 = 15.99491462  # 16O
m3 = 1.00782503   # proton
m4 = 16.9991315   # 17O
Elab = 20.0  # MeV

u = 931.494  # MeV/c^2 per amu

# CM kinematics
Ecm = Elab * m2 / (m1 + m2)
Q = (m1 + m2 - m3 - m4) * u  # Q-value in MeV
Ef_cm = Ecm + Q  # total CM KE in exit channel

# CM momenta squared (p^2 = 2*mu*E)
mu_i = m1 * m2 / (m1 + m2) * u  # MeV/c^2
mu_f = m3 * m4 / (m3 + m4) * u

ki2 = 2 * mu_i * Ecm
kf2 = 2 * mu_f * Ef_cm
ki = np.sqrt(ki2)
kf = np.sqrt(kf2)

# CM velocity (in momentum units: p_CM = m2*p_lab/(m1+m2) = m2*sqrt(2*m1*u*Elab)/(m1+m2))
# ... actually V_CM * m_total = p_lab in NR
# p_lab = sqrt(2 * m1 * u * Elab)
p_lab = np.sqrt(2 * m1 * u * Elab)
V_CM_p = p_lab / ((m1 + m2) * u)  # V_CM in c (natural units)
# ejectile CM speed: v3_CM = kf / (m3*u)
v3_CM = kf / (m3 * u)

gamma = V_CM_p / v3_CM  # dimensionless ratio
print(f"Ecm = {Ecm:.4f} MeV, Q = {Q:.4f} MeV, Ef_cm = {Ef_cm:.4f} MeV")
print(f"ki = {ki:.4f}, kf = {kf:.4f} (MeV/c)")
print(f"gamma = V_CM/v3_CM = {gamma:.6f}")

# --- Load C++ CM data ---
cpp_tm = np.loadtxt('/tmp/cpp_tmatch_dcs.txt')
cpp_wr = np.loadtxt('/tmp/cpp_wron_dcs.txt')

# Convert CM → lab for the outgoing proton
# theta_lab = arctan(sin(theta_CM) / (cos(theta_CM) + gamma))
# Jacobian: dOmega_CM/dOmega_lab = |d(cos theta_CM)/d(cos theta_lab)|
# = (1 + 2*gamma*cos(theta_CM) + gamma^2)^(3/2) / |1 + gamma*cos(theta_CM)|

theta_cm_rad = np.radians(cpp_tm[:, 0])
cos_cm = np.cos(theta_cm_rad)
sin_cm = np.sin(theta_cm_rad)

# Lab angle
theta_lab_rad = np.arctan2(sin_cm, cos_cm + gamma)
# Handle negative angles (arctan2 can give negative for back angles)
theta_lab_rad = np.where(theta_lab_rad < 0, theta_lab_rad + np.pi, theta_lab_rad)
theta_lab_deg = np.degrees(theta_lab_rad)

# Jacobian: dσ/dΩ_lab = dσ/dΩ_CM * (dΩ_CM/dΩ_lab)
# dΩ_CM/dΩ_lab = (1 + 2γ cos θ_CM + γ²)^(3/2) / |1 + γ cos θ_CM|
jacobian = (1 + 2*gamma*cos_cm + gamma**2)**1.5 / np.abs(1 + gamma*cos_cm)
dcs_lab_tm = cpp_tm[:, 1] * jacobian

# Same for Wronskian
theta_cm_rad_wr = np.radians(cpp_wr[:, 0])
cos_cm_wr = np.cos(theta_cm_rad_wr)
sin_cm_wr = np.sin(theta_cm_rad_wr)
theta_lab_rad_wr = np.arctan2(sin_cm_wr, cos_cm_wr + gamma)
theta_lab_rad_wr = np.where(theta_lab_rad_wr < 0, theta_lab_rad_wr + np.pi, theta_lab_rad_wr)
theta_lab_deg_wr = np.degrees(theta_lab_rad_wr)
jacobian_wr = (1 + 2*gamma*cos_cm_wr + gamma**2)**1.5 / np.abs(1 + gamma*cos_cm_wr)
dcs_lab_wr = cpp_wr[:, 1] * jacobian_wr

# --- Load Fortran lab data ---
ftn_dict = {}
with open('/tmp/ftn_full_output.txt') as f:
    in_section = False
    for line in f:
        if 'COMPUTATION OF CROSS SECTIONS' in line:
            in_section = True
            continue
        if in_section and 'TOTAL:' in line:
            in_section = False
        if in_section:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    angle = float(parts[0])
                    dcs = float(parts[1])
                    if 0 <= angle <= 180 and dcs > 0:
                        ftn_dict[angle] = dcs
                except ValueError:
                    pass
ftn_data = np.array(sorted(ftn_dict.items()))

print(f"\n0° check: C++ CM={cpp_tm[0,1]:.3f}, C++ lab={dcs_lab_tm[0]:.3f}, Ftn lab={ftn_data[0,1]:.3f}")
print(f"  ratio = {dcs_lab_tm[0]/ftn_data[0,1]:.4f}")
print(f"180° check: C++ CM={cpp_tm[-1,1]:.4f}, C++ lab={dcs_lab_tm[-1]:.4f}, Ftn lab={ftn_data[-1,1]:.4f}")
print(f"  ratio = {dcs_lab_tm[-1]/ftn_data[-1,1]:.4f}")

# --- Compute ratio at Fortran angles ---
ftn_angles = ftn_data[:, 0]
cpp_lab_interp_tm = np.interp(ftn_angles, theta_lab_deg, dcs_lab_tm)
cpp_lab_interp_wr = np.interp(ftn_angles, theta_lab_deg_wr, dcs_lab_wr)
ratio_tm = cpp_lab_interp_tm / ftn_data[:, 1]
ratio_wr = cpp_lab_interp_wr / ftn_data[:, 1]

pct_err_tm = (ratio_tm - 1) * 100
pct_err_wr = (ratio_wr - 1) * 100
print(f"\nTMATCH: mean |err| = {np.mean(np.abs(pct_err_tm)):.2f}%, max |err| = {np.max(np.abs(pct_err_tm)):.2f}%")
print(f"Wronskian: mean |err| = {np.mean(np.abs(pct_err_wr)):.2f}%, max |err| = {np.max(np.abs(pct_err_wr)):.2f}%")

# --- Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1]})

ax1.semilogy(theta_lab_deg, dcs_lab_tm, 'b-', linewidth=2, label='C++ TMATCH (CM→lab)')
ax1.semilogy(theta_lab_deg_wr, dcs_lab_wr, 'g--', linewidth=1.5, label='C++ Wronskian (CM→lab)')
ax1.semilogy(ftn_data[:, 0], ftn_data[:, 1], 'ro', markersize=5, label='Fortran Cleopatra (lab)')

ax1.set_xlabel('Lab Angle (deg)', fontsize=13)
ax1.set_ylabel('dσ/dΩ_lab (mb/sr)', fontsize=13)
ax1.set_title('¹⁶O(d,p)¹⁷O(5/2⁺ g.s.)  Elab = 20 MeV\nC++ CM→lab converted vs Fortran lab', fontsize=14)
ax1.legend(fontsize=12, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 180)

# Bottom: % error
ax2.plot(ftn_angles, pct_err_tm, 'b-o', markersize=3, linewidth=1.5, label=f'TMATCH (mean |err|={np.mean(np.abs(pct_err_tm)):.1f}%)')
ax2.plot(ftn_angles, pct_err_wr, 'g--s', markersize=3, linewidth=1.5, label=f'Wronskian (mean |err|={np.mean(np.abs(pct_err_wr)):.1f}%)')
ax2.axhline(y=0, color='r', linestyle='-', alpha=0.5)
ax2.fill_between(ftn_angles, -5, 5, alpha=0.1, color='green')
ax2.set_xlabel('Lab Angle (deg)', fontsize=13)
ax2.set_ylabel('% Error vs Fortran', fontsize=13)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 180)
ax2.set_ylim(-20, 20)

plt.tight_layout()
plt.savefig('/home/node/.openclaw/workspace/plots/o16dp_dcs_lab_comparison.png', dpi=150)
print("\nSaved to workspace/plots/o16dp_dcs_lab_comparison.png")
