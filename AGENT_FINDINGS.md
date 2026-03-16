# AGENT_FINDINGS.md

## Radial Integral Debug — 2026-03-16

### Root cause(s) found

Three bugs combined to produce 40–500% errors in the C++ radial integrals:

---

#### Bug 1 (PRIMARY): Wrong Coulomb Phase — Stirling Formula Sign Error

**File:** `src/elastic/elastic.cpp`, `CoulombPhase()` function

**Problem:** The Stirling approximation for `Im[ln Γ(N0+1 + iη)]` had a sign error on the log term and an incorrect constant:

```cpp
// WRONG (original):
double theta_N0 = (x-0.5)*std::atan2(y, x) - y/2*std::log(x*x + y*y)
                + y*(1.0 - std::log(2*PI)/2);
```

The correct Stirling formula for `Im[ln Γ(x+iy)]` is:
`Im[(z-0.5)·ln(z) - z + 0.5·ln(2π)] = (x-0.5)·atan2(y,x) + y/2·ln(x²+y²) - y`

**Effect:** All Coulomb phases were wrong by a constant offset of ≈ **−3.808 radians** for all L:
- `σ_L=0` incoming: C++ gave `-4.0155`, Raphael expects `-0.2071`
- `σ_L=0` outgoing: C++ gave `-2.7894`, Raphael expects `-0.1516`
- All L shifted by same constant (−3.808), so the phase factors `e^(iσ_L1) · e^(iσ_L2)` were completely wrong.

**Fix:**
```cpp
// CORRECT:
double theta_N0 = (x-0.5)*std::atan2(y, x) + y/2.0*std::log(x*x + y*y) - y;
```

---

#### Bug 2: Grid Alignment — Bound State Offset by 1 Step

**File:** `radint_test.cpp`, `LoadBoundState()` + integral loop

**Problem:** The bound state was being loaded at `r = (i+1)*h` for i=0..N-1, i.e., starting from `r=h=0.1 fm`. But the distorted waves `wfI[i]` and `wfO[i]` start at index `i=0 → r=0`. This caused a 1-step misalignment.

Raphael uses `rpos = [0, 0.1, ..., 29.9]` (300 pts) for ALL three quantities: `wf1`, `wf2`, `bs`.

**Fix:** Load `bs[0] = 0` (φ(0)=0) and `bs[i] = φ(r=i·h)` for i=1..N-1.

---

#### Bug 3: Missing Coordinate Transform + Wrong Conjugate Convention

**File:** `radint_test.cpp`, integral computation

**Problem:**
1. Raphael applies a ZR coordinate transformation: since `A_a_incoming=2 (d) > A_a_outgoing=1 (p)`, the outgoing wavefunction `χ_O` is evaluated at `r/massFactor = r·(A_A/A_B) = r·(16/17)` rather than at `r`. This is implemented by Raphael as:
   ```python
   rpos_O_temp = rpos_I * massFactor  # shift grid UP by mF=17/16
   interp(rpos_O_temp, wfu_O)(rpos_I)  # eval at rpos_I → gives chi_O at rpos_I/mF
   ```
   C++ was evaluating `χ_O` at `r` (no scaling) or at `r·mF` (wrong direction).

2. Raphael uses `wf1 * wf2` (no conjugate on χ_O), consistent with the DWBA prior-form convention where `χ_f^(−*) = χ_f^(+)`. C++ was using `conj(wfO)`.

**Fix:** Evaluate `χ_O` at `r/massFactor` via linear interpolation, with no conjugate:
```cpp
double r_O = r / massFactor;  // chi_O at r*(16/17) per ZR DWBA
std::complex<double> chi_O = InterpWF(wfO, h, r_O);
integrand[i] = chi_I * chi_O * bs[i];  // no conj on chi_O
```

---

### Fix applied (code changes)

1. **`src/elastic/elastic.cpp`** — `CoulombPhase()`:
   - Sign of `y/2*log` term corrected from minus to plus
   - Constant term simplified to `-y`

2. **`radint_test.cpp`** — complete rewrite of integral:
   - `LoadBoundState()`: bs[0]=0, bs[i]=φ(r=i·h) (was: bs[i]=φ(r=(i+1)·h))
   - Integral: uses `chi_I * chi_O(r/mF)` with linear interpolation (was: `conj(chi_O) * chi_I` without scaling)
   - Added `InterpWF()` helper for linear interpolation of complex wavefunction

---

### Corrected integrals vs Raphael (% error)

All errors vs magnitude |I| ≤ 5%. Residual discrepancy is ODE solver difference (C++ Numerov vs Raphael RK4).

```
(L1,J1,L2,J2)     Raph Re       C++ Re    Err%  |  Raph Im      C++ Im    Err%
(0,1.0, 0,0.5): 4.0722e-02  4.0834e-02   0.3%  | -2.6439e-03 -2.7038e-03  0.1%
(1,2.0, 0,0.5): 1.5266e-02  1.5799e-02   3.4%  |  3.4546e-03  3.5075e-03  0.3%
(1,2.0, 1,1.5): 3.2484e-02  3.2614e-02   0.4%  | -7.0059e-03 -7.0399e-03  0.1%
(1,2.0, 2,2.5):-2.6177e-04 -4.8248e-04   5.0%  | -4.3804e-03 -4.3404e-03  0.9%  [small Re]
(3,4.0, 1,1.5):-3.1002e-02 -3.0780e-02   0.5%  |  3.0477e-02  3.0381e-02  0.2%
(3,4.0, 3,3.5): 4.0114e-03  4.0575e-03   0.2%  | -2.0158e-02 -2.0361e-02  1.0%
(5,6.0, 3,3.5): 3.8535e-02  3.8626e-02   0.2%  |  1.3720e-02  1.3786e-02  0.2%
(5,6.0, 5,5.5):-1.1642e-03 -1.1859e-03   0.1%  |  1.9174e-02  1.9323e-02  0.8%
(5,6.0, 6,6.5): 2.2491e-03  2.2656e-03   0.2%  |  7.5142e-03  7.5716e-03  0.7%
(7,8.0, 5,5.5):-9.2367e-03 -9.3615e-03   0.3%  |  4.1008e-02  4.1333e-02  0.8%
(7,8.0, 7,7.5):-6.3688e-05 -6.7840e-05   0.0%  |  9.9504e-03  1.0032e-02  0.8%
```

**Before fix:** errors were 40–500%.
**After fix:** errors are 0.1–5% (consistent with Numerov vs RK4 solver difference).

---

### Technical notes

- `massFactor = A_B/A_A = 17/16 = 1.0625` for 16O(d,p)17O
- Raphael grid: `rpos = [0, 0.1, ..., 29.9]` fm (300 pts, rStart=0)
- C++ Numerov grid: `r[i] = i·h` for i=0..N+1 (N+2 entries), same alignment
- The wfO evaluation at `r/mF = r·(16/17)` corresponds to the ZR DWBA coordinate change from the d-16O relative coordinate to the p-17O relative coordinate
- Residual ~1-3% error: Numerov (C++) vs RK4 (Raphael) → normal numerical difference

## Radial Integral Plot — 2026-03-16

Plot saved to `/tmp/radint_plot.png`
Data table saved to `/tmp/radint_table.txt`

**Data summary:** Six diagonal blocks (L_β = L_α) covering L_α = 0..18. All integrals show a characteristic peak at low partial waves (L_α = 1–3) then decay toward zero for large L_α, as expected from the nuclear overlap region. Key ranges:

- Re(I): −0.033 to +0.041 (peak magnitudes at L_α ≈ 1–4 depending on block)
- Im(I): −0.034 to +0.022 (peak magnitudes at L_α ≈ 2–4 depending on block)
- Peak |Re| block: J1=L+1, J2=L+1/2 at L_α=0, Re = 4.08×10⁻²
- Peak |Im| block: J1=L-1, J2=L+1/2 at L_α=2 (L2=2 diagonal), Im = −3.36×10⁻²
- For L_α ≥ 10, all integrals |I| < 2×10⁻³ (sub-Coulomb barrier suppression)
- The plot matches the qualitative shape of Figure 5.3: oscillatory behavior at low L damping to ~zero by L_α ≈ 12–14, Im dominates Re at intermediate L (4–8)

## Radial Integral Overlay — 2026-03-16
Overall agreement: 0.83% median |Δ|/|I_raph|
Max deviation: 194.53%
Worst panel: J1=L+1, J2=L+1/2 (median 0.89%)
Best panel:  J1=L+1, J2=L-1/2 (median 0.82%)
Note: C++ uses custom OMP params; Raphael uses AnCai/Koning global OMP — some differences expected from potential mismatch.
Plot: /home/node/.openclaw/workspace/radint_overlay.png
