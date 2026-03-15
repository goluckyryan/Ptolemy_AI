# 60Ni(d,d) Elastic Scattering — Session Handoff

## Goal
Match C++ elastic DCS to Fortran Ptolemy for ⁶⁰Ni(d,d) at Elab=60 MeV (30 MeV/u), spin-1 deuteron.

## Current Status

### What works ✅
- **Proton case (spin-½):** `run_o16pp_cpp.cpp` → ¹⁶O(p,p) matches Ptolemy to <0.01% at ALL angles
- **Spin-½ S-matrix:** Perfect match (factor-4 + sign fix in `potential_eval.cpp`)
- **Spin-1 S-matrix at low L:** Matches Ptolemy at L=1-6 to ~1%

### What doesn't work ❌
- **Deuteron DCS:** C++ agrees ~5-10% at forward angles, diverges to ~60% at backward angles
- **Spin-1 S-matrix at large L:** C++ diverges from Ptolemy by 5-20% for L≥7

## Files
- **C++ source:** `/home/node/working/ptolemy_2019/Cpp_AI/Raphael_AI/run_ni60dd_cpp.cpp`
- **Potential eval:** `/home/node/working/ptolemy_2019/Cpp_AI/src/dwba/potential_eval.cpp`
- **Ptolemy output:** `ni60_dd_ptolemy.out` (reference)
- **Data:** `data_ni60dd_ptolemy.txt`, `data_ni60dd_cpp.txt`, `data_ni60dd_raphael.txt`

## Key Physics Parameters
- An & Cai (2006) global deuteron OMP at Ed=60 MeV for d+⁶⁰Ni:
  - V=77.917 r0=1.1500 a=0.7222
  - Wi=4.836  r0=1.3305 a=0.8295
  - Wsi=8.994 r0=1.3728 a=0.5468
  - Vso=3.557 r0=0.9720 a=1.011
  - RC=1.303
- Kinematics: Ecm≈58.06 MeV, k≈2.319 fm⁻¹, eta≈0.805, Lmax=44

## Current potential_eval.cpp SO formula
```cpp
V_so_real = -pot.VSO * 4.0 * (1.0/r) * deriv;  // negated, factor 4
// deriv = -(1/ASO) * exp/(1+exp)^2  (negative)
// So V_so_real = +4*VSO*exp/(ASO*r*(1+exp)^2)  → POSITIVE for VSO>0
```

## Current wavelj_loc Lambda formula
```cpp
Lambda = (J*(J+1) - L*(L+1) - S*(S+1)) / 2.0;  // standard L·S
// For S=1/2: Lambda = [-3/4 + J(J+1) - L(L+1)] / 2
// For S=1:   Lambda = [-2   + J(J+1) - L(L+1)] / 2
```

## Root Cause Hypothesis (confirmed from Fortran WAVSET source)

Ptolemy's spin-orbit in WAVSET:
```
V_SO = VSO * 4 * L·S * (1/r) * d/dr[WS]
L·S = [J(J+1) - S(S+1) - L(L+1)] / 2
```

But WAVELJ implementation splits it as:
```fortran
SDOTL = (J(J+1) - S(S+1) - L(L+1)) / (2*S)  ← divides by 2S not 2!
Vso_array = WOODSX(type=2) = 2*VSO/r * d/dr[WS]
effective = SDOTL * Vso_array
```

For spin-½ (S=½, 2S=1): SDOTL × 2/r = L·S/½ × 2/r = 4·L·S/r ✓
For spin-1  (S=1,  2S=2): SDOTL × 2/r = L·S/1  × 2/r = 2·L·S/r ≠ 4·L·S/r ✗

**Bottom line: Ptolemy actually uses `4*L·S/r*df/dr` for spin-½ but `2*L·S/r*df/dr` for spin-1 (due to the `/JSPS` division).**

## Things Already Tried (didn't fix it)
1. ❌ Factor 2 → 4 with negation (fixes spin-½, breaks spin-1 slightly)
2. ❌ Lambda = L·S/S (fixes spin-1 S-matrix magnitude pattern but breaks spin-½)
3. ❌ Matching radius (30 fm vs 50 fm): no effect
4. ❌ Finer step h=0.05: no effect

## Key Insight Still Needed
The S-matrix discrepancy at L≥7 cannot be just the SO formula — even with perfectly-matched SO, L≥7 has |S|≈1 (transparent). The real S-matrix spread at L=1-3 is only ~0.3% between J=L-1 and J=L+1. So even a large SO error shouldn't cause >60% DCS difference at backward angles.

**The DCS formula itself may be wrong for spin-1.** Ptolemy's elastic DCS for spin-1 is NOT just:
  dσ/dΩ = (1/3) Σ_{v0,v} |fc δ_{vv0} + fN(v,v0)|²

Ptolemy uses the LX-basis (LX=0,1,2 for spin-1) and computes amplitudes F(LX,MX), then:
  dσ/dΩ = Σ_{LX,MX} FMNEG × |F(LX,MX)|² × 10

where FMNEG=1 for MX=0, 2 for MX>0. This is equivalent to the (v,v0) sum, BUT only if the Coulomb amplitude is added correctly to EACH separate (LX,MX) component — not to the full amplitude.

## What To Try Next

### Option A: Fix the DCS formula
Check whether Raphael's `DCSUnpolarized` is actually correct for spin-1.
Compare term-by-term at θ=40° between Raphael and C++:
- Print G(v,v0,L) for each L, v0, v
- Print nuclear amplitude fN(v,v0) at θ=40°
- Print Coulomb fc
- Print |fc + fN|² vs Ptolemy's |F_LX0|² + 2|F_LX1|² + 2|F_LX2|²

### Option B: Use S-matrix from Ptolemy directly in C++ DCS formula
We have Ptolemy's LX-basis S-matrix from ni60_dd_ptolemy.out.
Convert to J-basis using T_inv (validated correct).
Apply C++ DCS formula with those S values.
If DCS matches Ptolemy → DCS formula correct, S-matrix wrong.
If DCS still wrong → DCS formula wrong.

### Option C: Spin-0 sanity check
Run Ptolemy with VSO=0, run C++ with twoS=0.
Should get spinless DCS — compare. This isolates SO issues from formula issues.

## T-matrix (LX → J basis) for spin-1, validated
```python
T_inv = np.array([
    [1.0,        -1.73205081,  2.23606798],   # J=L-1
    [1.0,        -0.8660254,  -1.11803399],   # J=L
    [1.0,         0.8660254,   0.2236068 ]    # J=L+1
])
# SJ = T_inv @ [SLX0, SLX1, SLX2]
```

## Numerical Reference (from Ptolemy, ni60_dd_ptolemy.out)
Ptolemy DCS (σ/σ_Ruth) at key angles:
  20°: 0.7004,  40°: 0.4071,  60°: 1.654,  80°: 0.7113
 100°: 0.2108, 120°: 0.0805, 140°: 0.0437, 160°: 0.0345
 175°: 0.0868, 180°: 0.0728

Ptolemy S-matrix (LX-basis) for first few L:
L=0: LX0 = -0.23018 - 0.00449i
L=1: LX0 = -0.034704+0.22597i  LX1=-0.019032-0.002117i  LX2=3.26e-5-4.46e-4i
L=2: LX0 =  0.15207 +0.16814i  LX1=-0.022615+0.022908i  LX2=-0.001119-0.001024i
L=3: LX0 =  0.22018 +0.02033i  LX1=-0.002156+0.044891i  LX2=-0.003098-3.84e-5i

## RAG Tools (use before any grep/read)
```
python3 /home/node/working/ptolemy_2019/rag/search.py "query" 2>/dev/null
python3 /home/node/working/ptolemy_2019/rag_cpp/search_cpp.py "query" 2>/dev/null
```

## Key Fortran source locations
- WAVSET SO definition: source.mor line ~31960 (WAVSET comments)
- WAVELJ SO application: source.mor line ~30668 (SDOTL = ...)
- WOODSX type-2 formula: source.mor line ~32610 (SUBROUTINE WOODSX)
- ELDCS elastic DCS: source.mor line ~12989 (SUBROUTINE ELDCS)
- BETCAL amplitude: source.mor line ~3358 (SUBROUTINE BETCAL)
- AMPCAL: source.mor line ~220 (SUBROUTINE AMPCAL)

---

## Raphael solveSE.py Bug Fixes (2026-03-15 evening) — CONFIRMED

### Bugs fixed (Ryan approved)
1. **L·S formula**: `S*(S)` → `S*(S+1)` in `solveSE.py` LS() method
2. **SO scale**: Added `so_scale = 1/(2*S)` to SO term in `__PotentialValue()`

### Verification
- 60Ni(d,d) S-matrix: ALL LX=0,1,2 match Ptolemy <0.5% for L=0..15 ✅
- 60Ni(d,d) DCS: <0.3% vs Ptolemy at all angles 5°-180° ✅

### ⚠️ WARNING: spin-½ impact needs re-check
Fix 1 changes LS() for S=½ by −0.125 per partial wave (S²=0.25 vs S(S+1)=0.75).
Must re-verify 16O(p,p) and 60Ni(p,p) still match Ptolemy after this change.
