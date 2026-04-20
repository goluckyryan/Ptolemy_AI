# Ptolemy Zero-Range DWBA — Architecture & Physics

*Written 2026-04-19 after comparison of Cleopatra NZRDWBA vs Raphael ZR*

---

## 1. What NZRDWBA Actually Does

Despite the name "zero-range", Ptolemy's `NZRDWBA` does **not** use a true Dirac delta function for the neutron-proton interaction. Instead:

```
VPHI(r_I) = φ_AV18(r_np(r_I)) × V_np(r_np(r_I))
```

where `φ_AV18` is the actual AV18 deuteron wavefunction evaluated at each point on the 2D integration grid. The integral is still done over (RI, RO) — a full 2D grid with ~1600 points.

The "zero-range" refers to the use of the **D0 normalization constant** (D0 = 1.55×10⁴ MeV²·fm³) rather than computing the full finite-range overlap. The AV18 wavefunction provides the spatial shape; D0 provides the absolute normalization at r_np=0.

**Consequence:** Cleopatra NZRDWBA ≠ true ZR (Raphael/textbook). Cleopatra is better described as "ZR-normalized finite-range with AV18 shape".

---

## 2. Coordinate System

For A(a,b)B transfer (e.g. 206Hg(d,p)207Hg):

### Jacobi Coordinates
- `r_I` = incoming Jacobi vector (deuteron CM relative to target A)
- `r_O` = outgoing Jacobi vector (proton relative to residual B)
- `r_np` = neutron-proton separation

### ZR Constraint
In zero-range approximation, `r_np = 0`, which in Ptolemy's (RI, RO) coordinates gives:

```
S2 × r_I + T2 × r_O = 0
→  r_O = (A_A / A_B) × r_I
```

For 206Hg(d,p)207Hg (A_A=206, A_B=207):
```
r_O = (206/207) × r_I = 0.9952 × r_I
```

The outgoing coordinate is **slightly compressed** relative to incoming (always < 1 for stripping).

### GRDSET Coordinate Coefficients
Computed from mass ratios `BRATMS(1) = AMX/AMB`, `BRATMS(2) = AMX/AMBIGA`:

```fortran
TEMP = 1 / (BRATMS(1) + BRATMS(2)*(1+BRATMS(1)))
S1 = (1+BRATMS(1)) * (1+BRATMS(2)) * TEMP  ! ≈ 104 for 206Hg(d,p)
T1 = -(1+BRATMS(2)) * TEMP                  ! ≈ -103.5
S2 = (1+BRATMS(1)) * TEMP                   ! ≈ 103.5
T2 = -S1                                    ! ≈ -104
JACOB = S1³                                 ! Jacobian of coord transform
```

ZR line: `r_O = -(S2/T2) × r_I = (S2/S1) × r_I = A_A/A_B × r_I`

---

## 3. Mass Rescaling in Radial Integral

For the 1D ZR integral (along the ZR constraint line), the integration variable is `r_I`. The outgoing distorted wave must be evaluated at `r_O = (A_A/A_B) × r_I`.

**Correct massFactor for stripping (A_B > A_A):**
```
massFactor = A_A / A_B  (< 1 for stripping)
```

**Correct massFactor for pickup (A_A > A_B):**
```
massFactor = A_B / A_A  (< 1 for pickup)
```

In both cases, `massFactor = min(A_A,A_B) / max(A_A,A_B) < 1`.

The Jacobian factor `dr_O = massFactor × dr_I` enters the radial integral as an extra `massFactor` multiplicative factor.

**Bug fixed in Raphael (2026-04-19):** Original code used `A_B/A_A` for stripping (> 1, wrong). Fixed to `A_A/A_B` (< 1, correct).

---

## 4. Spin-Orbit Convention: Ptolemy vs Raphael

Ptolemy computes the LS expectation value as:

```fortran
SDOTL = (J(J+1) - L(L+1) - S(S+1)/4) / (2S)   ! divided by 2S internally
```

So the SO Hamiltonian applied is:
```
H_SO = V_SO_input × f(r) × (L·S) / (2S)
```

For the **physical** SO Hamiltonian `H = V_SO_phys × f(r) × (L·S)`, this means:
```
V_SO_phys = V_SO_ptolemy / (2S)
```

Raphael uses the physical convention (no 2S division), so to match Ptolemy potentials:
- **Deuteron (S=1):** `V_SO_raphael = V_SO_ptolemy / 2`
- **Proton (S=1/2):** `V_SO_raphael = V_SO_ptolemy` (no change, 2S=1)
- **Bound neutron (S=1/2):** no correction needed

*Note: in practice for 206Hg(d,p), this 2S correction has negligible effect on the DCS shape (~0.1% change in mean ratio).*

---

## 5. Bound State Parameters

Ptolemy's TARGET block uses `r0target` convention:
```
R0 = r0 × A_target^{1/3}   (target nucleus only, not projectile+target)
```

For 207Hg bound state (1g9/2, nodes=1, l=4, jp=9/2):
```
r0=1.25 fm,  a=0.65 fm  → R0 = 1.25 × 207^{1/3} = 7.382 fm
vso=6.0 MeV, rso0=1.10 fm → RSO = 1.10 × 207^{1/3} = 6.497 fm
rc0=1.30 fm
BE = -3.3445 MeV  (from Ptolemy Q-value, not AME2020 = -3.613 MeV)
V_central (fitted) = -45.654 MeV
KAPPA = 0.40077 fm^-1
```

**Node counting**: Ptolemy `nodes=1` = Raphael node=1 (one sign change of u(r) within r_cut). For 1g9/2 with Ptolemy BE, V_central is found at ~-45.5 to -45.7 MeV.

---

## 6. D0 Normalization

For (d,p) stripping:
```
D0 = 1.55 × 10^4  MeV²·fm³
A_lsj² = D0 × (2s_d+1)/(2s_n+1) = D0 × 3/2 = 2.325 × 10^4
```

D0 is a purely multiplicative constant — it affects overall DCS magnitude but **not** the angular shape.

---

## 7. Raphael vs Cleopatra ZR Comparison (206Hg(d,p)207Hg)

With matched kinematics (same E_cm_out, k_out, mu_out, BE, bound state potential):

| Quantity | Agreement |
|---|---|
| Elastic S-matrix |S| | ~1% (L=0: 0.01%, L=1: 8%, L>3: <1%) |
| Bound state V_central | 0.26% (45.53 vs 45.65 MeV) |
| DCS mean ratio (Raphael/Cleopatra) | 1.013 ± 0.066 |
| DCS at 0° | 1.28× (Raphael larger) |
| DCS at 90° | ~0.99 |
| DCS at 180° | ~0.94 |

**Root cause of shape difference**: Ptolemy integrates with actual AV18 φ(r_np) over 2D grid; Raphael uses true D0·δ(r_np). The AV18 finite width (~2 fm) enhances forward angles in Ptolemy, causing Raphael to be 28% larger at 0°.

---

## 8. Angular Momentum Coupling in ZR Beta Function

In Raphael's `AngDist`, the summation index `m` is:
```python
for m in np.arange(-j + mb - ma, j + mb - ma + 1, 1)
```

For (d,p): spin_a=1 (integer), spin_b=1/2 (half-integer) → `mb-ma` is half-integer → `m` takes **half-integer values**. These are unphysical (L2 is integer, requires integer m for Y_L^m). However, they produce **zero** automatically because `PreCalClebschGordan` only stores CG values for physically valid quantum numbers — half-integer m for integer L2 is never computed → array lookup returns 0.

**Bottom line**: Raphael's ZR result is numerically correct despite the half-integer m loop entries. The angular algebra is sound.

---

## 9. Key Findings Summary

1. ✅ **massFactor fix**: `A_A/A_B` (not `A_B/A_A`) for stripping — fixes Raphael's ZR coordinate
2. ✅ **Ptolemy NZRDWBA is NOT true ZR**: uses AV18 φ(r_np) shape over 2D grid
3. ✅ **2S factor**: Ptolemy divides VSO by 2S internally; negligible effect on DCS for this reaction
4. ✅ **Binding energy matters**: AME2020 (-3.613 MeV) vs Ptolemy Q-value (-3.3445 MeV) shifts DCS shape
5. ✅ **D0 is purely normalization**: doesn't affect angular shape
6. ✅ **S-matrix matches**: Raphael and Cleopatra elastic S-matrices agree to <2% for L>1
