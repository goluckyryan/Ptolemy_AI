#!/usr/bin/env python3
"""
Extract S-matrix from Raphael for 60Ni(d,d) at 60 MeV, compare with Ptolemy LX=0,1,2.
Uses JPTOLX formula from source.mor. Ptolemy output parsed with correct sign handling.
"""
import sys, re
import numpy as np
sys.path.insert(0, '/home/node/working/PtolemyGUI/Raphael')

from distortedWave import DistortedWave
from solveSE import WoodsSaxonPot, WS_SurfacePot, SpinOrbit_Pot, CoulombPotential
from clebschGordan import sixj

# ── Parse Ptolemy output ────────────────────────────────────────────
ptol = {}
with open('/home/node/working/ptolemy_2019/Cpp_AI/Raphael_AI/ni60_dd_ryan_ptolemy.out') as f:
    for line in f:
        m = re.match(r'\s+(\d+)\s+(\d+)\s+(\d+)\s+([+-]?[\d.E+-]+)\s+(\+-?|\+)\s*([\d.E+-]+)\s+I', line)
        if m:
            L, LX = int(m.group(1)), int(m.group(3))
            Re = float(m.group(4))
            Im = float(m.group(6))
            if m.group(5) == '+-': Im = -Im
            ptol[(L, LX)] = complex(Re, Im)

# ── Raphael: Ryan's exact params ────────────────────────────────────
kaka = DistortedWave("60Ni", "d", 60)
kaka.ClearPotential()
kaka.AddPotential(WoodsSaxonPot(-81.919,       1.150, 0.768), False)
kaka.AddPotential(WoodsSaxonPot(-4.836j,       1.330, 0.464), False)
kaka.AddPotential(WS_SurfacePot(-8.994j,       1.373, 0.774), False)
kaka.AddPotential(SpinOrbit_Pot(-3.557+0.000j, 0.972, 1.011), False)
kaka.AddPotential(CoulombPotential(1.303),                    False)
kaka.CalScatteringMatrix(True, verbose=False)

print("="*95)
print("Raphael vs Ptolemy — 60Ni(d,d) Elab=60 MeV (Ryan's exact An&Cai params)")
print("JPTOLX: S_LX = Σ_J (2J+1)/sqrt((2S+1)(2L+1)) * phase * sqrt(2LX+1) * W(L,L,S,S;LX,J) * S_J")
print("="*95)
S_spin = kaka.S   # =1 (deuteron)
JSPS   = int(2*S_spin)  # =2
print(f"  S={S_spin}, JSPS={JSPS},  Ecm={kaka.Ecm:.4f} MeV,  k={kaka.k:.5f} fm^-1,  eta={kaka.eta:.5f}")
print(f"  Ptolemy parsed: {len(ptol)} (L,LX) entries")
print()

sm = kaka.ScatMatrix

print(f"{'L':>3} {'LX':>3}  {'Raph Re':>12} {'Raph Im':>12}  {'|Raph|':>8}   {'Ptol Re':>12} {'Ptol Im':>12}  {'|Ptol|':>8}  {'Δ%':>7}")
print("-"*105)

for L in range(min(len(sm), 16)):
    L2 = 2*L
    # Build 2J -> S_J mapping (Raphael idx=0->J=L-S, ..., idx=2S->J=L+S)
    S_J = {}
    for idx in range(len(sm[L])):
        J = L - S_spin + idx
        if J >= 0:
            S_J[int(round(2*J))] = sm[L][idx]

    for LX in range(min(JSPS+1, 2*L+2)):
        LX2 = 2*LX
        S_LX = complex(0, 0)
        for JP2, Sv in S_J.items():
            # JPTOLX coefficients (Fortran source, all args are 2J integers)
            I    = (JSPS - JP2 + L + L) // 2
            sign = (-1)**I
            coef = sign * (JP2+1) / np.sqrt((JSPS+1)*(2*L+1))
            # RACAH(2L,2L,JSPS,JSPS,LX2,JP2) = (-1)^{(2L+2L+JSPS+JSPS)/2} * sixj(L,L,LX,S,S,J)
            ph_r = (-1)**((2*L + 2*L + JSPS + JSPS)//2)
            W    = ph_r * sixj(L, L, LX, JSPS/2, JSPS/2, JP2/2)
            S_LX += coef * np.sqrt(LX2+1) * W * Sv

        S_ptol = ptol.get((L, LX))
        if S_ptol is not None:
            pct  = abs(S_LX - S_ptol) / max(abs(S_ptol), 1e-6) * 100
            flag = "✅" if pct < 2 else ("⚠️" if pct < 10 else "❌")
            print(f"  {L:2d}  LX={LX}  {S_LX.real:+12.6f} {S_LX.imag:+12.6f}  {abs(S_LX):8.5f}   "
                  f"{S_ptol.real:+12.6f} {S_ptol.imag:+12.6f}  {abs(S_ptol):8.5f}  {pct:7.2f}% {flag}")
        else:
            print(f"  {L:2d}  LX={LX}  {S_LX.real:+12.6f} {S_LX.imag:+12.6f}  {abs(S_LX):8.5f}   (no Ptolemy ref)")

    # Individual J values
    for JP2, Sv in sorted(S_J.items()):
        print(f"       J={JP2//2}  {Sv.real:+12.6f} {Sv.imag:+12.6f}  {abs(Sv):8.5f}")
    print()

# ── Verify DCS still matches after fixes ─────────────────────────────
if __name__ == '__main__':
    print("\nVerifying DCS vs Ptolemy after S-matrix fix...")
    kaka.CalAngDistribution(180, 1.0, None, False)
    # Sample a few angles
    sample_angles = [5, 10, 20, 40, 60, 90, 120, 150, 180]
    for th in sample_angles:
        dcs = kaka.DCSUnpolarized(th)
        print(f"  theta={th:3d} deg: dσ/dΩ = {dcs:.4f} mb/sr")

# Quick DCS spot check
import importlib
kaka2 = DistortedWave("60Ni", "d", 60)
kaka2.ClearPotential()
kaka2.AddPotential(WoodsSaxonPot(-81.919,       1.150, 0.768), False)
kaka2.AddPotential(WoodsSaxonPot(-4.836j,       1.330, 0.464), False)
kaka2.AddPotential(WS_SurfacePot(-8.994j,       1.373, 0.774), False)
kaka2.AddPotential(SpinOrbit_Pot(-3.557+0.000j, 0.972, 1.011), False)
kaka2.AddPotential(CoulombPotential(1.303),                    False)
kaka2.CalScatteringMatrix(True, verbose=False)
kaka2.CalAngDistribution(180, 5.0, None, False)
print("\nDCS spot check (after both fixes):")
for th in [5,10,20,40,60,90,120,150,180]:
    dcs = kaka2.DCSUnpolarized(th)
    print(f"  theta={th:3d}:  {dcs:.4f} mb/sr")
