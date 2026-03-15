# Ptolemy DWBA: How It Works (and How It Compares to Your Handbook)

**Dudu's guide for Ryan 🐱**  
*Cross-referenced with: `handbook_of_direct_nuclear_reaction_for_retarded_theorist_v1-1.pdf` and `source.mor`*

---

## Overview: The Big Picture

Ptolemy computes the **finite-range DWBA differential cross section** for a transfer reaction A(a,b)B.
For 33Si(d,p)34Si: projectile a=d, ejectile b=p, transferred particle x=n, target A=33Si, residual B=34Si.

The calculation flows through these stages:

```
INPUT PARSING
    │
    ▼
BOUND STATES  ← solve Schrödinger eq. for n inside 34Si and inside d
    │
    ▼
DISTORTED WAVES  ← solve optical model Schrödinger eq. for d+33Si and p+34Si
    │
    ▼
COORDINATE SETUP (GRDSET)  ← build 2D grid in (r_α, r_β)
    │
    ▼
RADIAL INTEGRALS (INELDC)  ← integrate χ*(r_β) · K(r_α,r_β) · χ(r_α) d r_α d r_β
    │
    ▼
S-MATRIX ELEMENTS (SFROMI)  ← apply kinematic factor, CG coefficients, phase
    │
    ▼
CROSS SECTION (XSECTN+BETCAL)  ← coherent sum → dσ/dΩ
```

---

## Step 1: Bound States (`CalculateBoundState` / Ptolemy's WSAXON)

### What it does
Solve the radial Schrödinger equation for the transferred neutron:

$$\left[\frac{d^2}{dr^2} - \frac{l(l+1)}{r^2} + \frac{2\mu}{\hbar^2}\left(E - V_{WS}(r) - V_{SO}(r)\right)\right] u_{lj}(r) = 0$$

with **negative energy** E = −(binding energy), outward boundary condition `u(0) = 0`, and an exponentially decaying tail at large r.

**Two bound states are needed:**
- `φ_T(r)` = neutron bound in **34Si** (the target): orbital 0d3/2, binding energy ~7.55 MeV, L=2, j=3/2
- `φ_P(r)` = neutron bound in **deuteron** (the projectile): s-wave, binding energy 2.22 MeV, L=0

**Normalization:**
$$\int_0^\infty u_{lj}^2(r)\, dr = 1 \quad\Rightarrow\quad \phi_{lj}(r) = u_{lj}(r)/r$$

### How the depth V is found
Binary search / iterative refinement to find V such that the wave function has exactly `n` nodes and decays exponentially. Uses **outward + inward Numerov matching** at a classical turning point.

### Numerov method (your Handbook §8.7)
Your handbook derives the Runge-Kutta method. Ptolemy uses the related **Numerov method**, optimized for second-order ODEs with no first-derivative term:

$$u_{i+1} = \frac{2(1 - \frac{5h^2}{12}f_i)u_i - (1 + \frac{h^2}{12}f_{i-1})u_{i-1}}{1 + \frac{h^2}{12}f_{i+1}}$$

where $f_i = k^2 - \frac{l(l+1)}{r_i^2} - \frac{2\mu}{\hbar^2}V(r_i)$.

---

## Step 2: Distorted Waves (`WavElj` / Ptolemy's `WAVELJ`)

### What it does
Solve the **optical model** Schrödinger equation for the elastic channels:

$$\left[\frac{d^2}{dr^2} - \frac{l(l+1)}{r^2} + k^2 - \frac{2\mu}{\hbar^2}\left(V_{opt}(r) + V_{so}(r) + V_C(r)\right)\right] \chi_{LJ}(r) = 0$$

The potential `V_opt` is complex (Wood-Saxon + imaginary surface/volume), hence `χ` is complex.

**Boundary conditions (your Handbook eq. after Eq. 27):**
$$\chi_{LJ}(r) \xrightarrow{r\to\infty} \frac{e^{i\sigma_L}}{2i}\left[S_{LJ}H_L^+(kr) - H_L^-(kr)\right]$$

where `σ_L` is the Coulomb phase shift, `S_LJ` is the S-matrix element, and `H_L±` are Coulomb Hankel functions.

**Normalization:** The distorted wave is normalized so that it asymptotically matches Coulomb functions `F_L` and `G_L`. This is what `alpha = chi_last / u_N` achieves in WAVELJ.

### Why L>0 is tricky (the centrifugal barrier problem)
For large L, the centrifugal term `L(L+1)/r²` dominates at small r, causing the Numerov solution to grow exponentially:  
- At r = 0.1 fm, L=2: centrifugal ~ 600 fm⁻², giving growth factor ~e^(√600 × 0.1) ≈ 11 per step  
- Over ~36 steps, growth reaches 10^36 → **floating point overflow**

**Ptolemy's WAVELJ fix** (lines ~30470–30489 in source.mor):  
When `|u| > BIGNUM` (= 10^30), rescale all stored values by `1/|u|` **AND zero out the inner region** where `|u * 1/|u_max|| < STEPI`. This cleanly separates the exponentially-growing inner region (which is unphysical — the wave function should be zero there for the regular solution) from the physically-meaningful outer region.

### S-matrix extraction (Wronskian matching)
At the outer boundary, match the Numerov solution to a linear combination of `F_L + G_L`:
$$u_N = \alpha_R F_L(r_N) + \alpha_I G_L(r_N) + S_L(\alpha_I F_L(r_N) - \alpha_R G_L(r_N))$$

The two-point Wronskian system gives S_L directly. This is exact — no approximation.

---

## Step 3: Coordinate Setup (`GRDSET`)

### The coordinate transformation (crucial for finite-range!)
In a transfer reaction A(a,b)B, the two relevant coordinates are:
- **r_α** = relative coordinate d ↔ 33Si (entrance channel)
- **r_β** = relative coordinate p ↔ 34Si (exit channel)
- **r** = internal coordinate n ↔ p (inside the deuteron)
- **R** = n ↔ 33Si (inside 34Si)

These are related by mass-weighted combinations (your Handbook Fig. 9):

$$\vec{r}_\alpha = S_1\vec{R} + T_1\vec{r} \qquad \vec{r}_\beta = S_2\vec{R} + T_2\vec{r}$$

For **stripping (ISTRIP=1)**, from Ptolemy source.mor lines 15876–15882:

| Coefficient | Formula | For d+33Si→p+34Si (xA=33, xb=1) |
|---|---|---|
| S1 | (1+xA)(1+xb) / denom | ≈ 1.029 |
| T1 | −(1+xb) / denom | ≈ −0.030 |
| S2 | (1+xA) / denom | ≈ 0.971 |
| T2 | −S1 | ≈ −1.029 |

where `denom = xA + xb + xA*xb + 1`.

### The 2D integration grid
Ptolemy integrates over `(r_α, r_β)` on a 2D Cartesian grid. For each pair (r_α, r_β), it evaluates:
- The x-coordinate: `x = √(S1²r_α² + S2²r_β² + 2S1S2r_α r_β cos θ_x)` — distance n↔33Si
- The bound state value `φ_T(x)` uses the target potential at that x
- Similarly for the projectile coordinate

The phi angle `φ_x` is integrated separately (Gauss quadrature over the azimuthal angle).

---

## Step 4: Radial Integrals (`INELDC`)

### The finite-range T-matrix kernel

The radial integral for a given (L_α, J_α, L_β, J_β) combination is:

$$I_{L_\alpha J_\alpha L_\beta J_\beta}^{lsj} = \frac{4\pi}{k_\alpha k_\beta} \int_0^\infty \int_0^\infty \chi_{L_\beta J_\beta}^{(-)*}(k_\beta r_\beta)\, K_{lx}(r_\alpha, r_\beta)\, \chi_{L_\alpha J_\alpha}^{(+)}(k_\alpha r_\alpha)\, dr_\alpha\, dr_\beta$$

The **kernel** `K` is formed from the nuclear form factors:

$$K_{L_x}(r_\alpha, r_\beta) = \sum_{\phi_x} \phi_T^*(x)\, V_{bx}(r)\, \phi_P(r)\, P_{L_x}(\cos\theta_x)\, \sin\theta_x\, d\theta_x$$

This is the projection of the product of bound states × interaction onto angular momentum transfer `L_x`, integrated over the internal angle.

**Measure:** `dr_α dr_β` (no extra r² factors; those are already inside φ = u/r → the r factors cancel)

### INELDC loop structure
```
for LI (incoming L, even and odd parity passes):
  compute distorted waves χ_Li for incoming channel
  for LO (outgoing L):
    compute distorted waves χ_Lo for outgoing channel
    for Lx (angular momentum transfer):
      compute CG angular coupling coefficients (A12 routine)
      accumulate radial integral I(Li, Lo, Lx)
call SFROMI to convert I → S-matrix element
```

---

## Step 5: S-Matrix Elements (`SFROMI`)

This is the **most important normalization step** — and where our C++ code was missing a factor.

### The kinematic factor
```fortran
FACTOR = 2 * DSQRT( AKI * AKO / (ES(1) * ES(2)) )
```
- `AKI` = k_α (fm⁻¹), `AKO` = k_β (fm⁻¹)
- `ES(1)` = E_cm of incoming channel (MeV), `ES(2)` = E_cm of outgoing channel (MeV)
- For our reaction: FACTOR = 2√(1.311 × 1.07 / (18.8 × 24.2)) ≈ 0.074

### The conversion formula
```fortran
TEMP = FACTOR * ATERM(LXP+1) / DSQRT(2*LASI+1.D0)
ITEST = LASI + LASO + 2*LXP + 1
IF ( MOD(ITEST, 4) .GE. 2 ) TEMP = -TEMP
S_matrix_element += (TEMP + i*TEMP) * I(Li, Lo, Lx)
```

- `ATERM` = spectroscopic coefficient (CG + SF)  
- `1/√(2La+1)` = normalization for the angular momentum of distorted wave channel  
- Phase `i^(La + Lo + 2Lx + 1)` (from partial wave expansion conventions)

---

## Step 6: Cross Section (`XSECTN` / `BETCAL`)

### From S-matrix to β_Lx
Each S-matrix element is sorted by angular momentum transfer Lx.  
The beta coefficient for Lx is the coherent sum:

$$\beta_{L_x} = \sum_{L_i, L_o} S(L_i, L_o, L_x) \cdot C(L_i, 0, L_o, 0; L_x, 0)$$

### The differential cross section
From your Handbook eq. after Eq. 36, the final formula is:

$$\frac{d\sigma}{d\Omega} = \frac{\mu_\alpha \mu_\beta}{\hbar^4} \cdot \frac{k_\beta}{k_\alpha} \cdot \frac{2J_B+1}{(2J_A+1)(2s_a+1)} \cdot \left|\sum_{L_x} \beta_{L_x} P_{L_x}(\cos\theta)\right|^2 \cdot \text{[units: mb/sr]}$$

Prefactor in practical units (converting fm² → mb, factor of 10):

$$\text{Prefactor} = \frac{(\mu_\alpha \mu_\beta)^2}{4\pi^2 \hbar^4} \cdot \frac{k_\beta}{k_\alpha k_\alpha k_\beta} \times 10 \quad \text{[mb/sr per |β|²]}$$

---

## Comparison: Ptolemy vs. Your Handbook

| Aspect | Your Handbook (§5, §5.8) | Ptolemy Implementation |
|---|---|---|
| **Approximation** | Both zero-range AND finite-range derived | **Finite-range only** (that's the whole point) |
| **T-matrix form** | Eq. 33: `<χ_β Φ_B \| W_β \| χ_α Φ_A>` | Same. `W_β = V_bx − (U_β − U_α)` (prior form) |
| **Kernel expansion** | Eq. 34: `I_βα(r_α, r_β)` via spherical harmonics | Same. GRDSET builds this via Gaussian quadrature in φ_x |
| **Coordinate transforms** | Fig. 9 with S1, T1, S2, T2 | GRDSET lines 15876–15882, same formulas |
| **Bound states** | §5.4: WS potential, match to Whittaker fn | WSAXON: bisection + inward-outward Numerov matching |
| **Distorted waves** | Eq. after 27: `χ_LJ → (e^iσ/2i)(S·H⁺ − H⁻)` | WAVELJ: outward Numerov + Wronskian S-matrix extraction |
| **Radial integral** | Eq. 5.8: `∫∫ χ_β* F_lsj(r_β,r_α) χ_α dr_α dr_β` | INELDC: same, 2D Gauss quadrature |
| **Kinematic factor** | Your Eq. 36: `μ_αμ_β/ℏ⁴ · k_β/k_α` | SFROMI: `2√(k_α k_β / E_cm_α E_cm_β)` (equivalent) |
| **Angular coupling** | 9j-symbol + CG coefficients (Eqs. 34–36) | A12 subroutine: same 9j + CG structure |
| **Phase convention** | `i^(L_α − L_β − l)` in Γ factor (your Eq.) | `i^(La + Lo + 2Lx + 1)` in SFROMI (equivalent with different bookkeeping) |
| **Cross section** | Eq. 36: sum over j, m with \|β\|² | BETCAL + XSECTN: same, organized by L_x |
| **Normalization of χ** | χ → (e^iσ/2i)(S·H⁺ − H⁻) normalized | ALPHAR/ALPHAI in WAVELJ: same asymptotic normalization |
| **Spin-orbit in waves** | §2.9.2, Eq. 18: standard `Λ(L,s,J)` | Same, with `Λ = ½[J(J+1) − L(L+1) − s(s+1)]` |
| **Units** | fm² → multiply by 10 for mb | Same: `×10` in XSECTN prefactor |

### Key Differences

**1. Finite-range vs. zero-range**  
Your handbook demonstrates both; the differences between Ptolemy and DWUCK4 in your Fig. 11 are **entirely due to this**. Ptolemy does the full 2D radial integral `∫∫ dr_α dr_β`, while DWUCK4 uses D₀²δ(r) to collapse it to a 1D integral.

**2. Satchler's spectroscopic coefficient convention**  
Your handbook §5.2.1 notes:
- Satchler: `|A_lsj(Satchler)|² = D₀² (2s_a + 1) S_lsj^{JA JB}`
- Ptolemy: `A_lsj = D₀ √(S_lsj^{JA JB})` (absorbs the √(2s_a+1) differently)

**3. The FACTOR in SFROMI**  
Your handbook writes the kinematic factor as `μ_αμ_β/(ℏ²)²`. Ptolemy factors this as `2√(k_αk_β / E_cm_α E_cm_β)` times a separate prefactor in XSECTN. These are the same thing after expansion, but need to match exactly in the C++ code.

**4. Phase of the radial integral**  
Your handbook eq. for `Γ` has phase `i^(L_α − L_β − l)`.  
Ptolemy's SFROMI uses `i^(La + Lo + 2*Lx + 1)`.  
For our case (La=Lo=Lx=2): handbook → `i^0 = 1`, Ptolemy → `i^9 = i`. This **imaginary unit phase** must match exactly.

**5. S-matrix output format (your Appendix 8.3)**  
You noted: *"The scattering matrix output of Ptolemy for spin-1/2 is actually need additional manipulation... S_l± ≈ R₀ ± R₁."*  
This is correct — Ptolemy outputs the S-matrix **before the 9j-symbol recoupling**. The actual physical S_LJ is recovered via that recoupling.

---

## C++ Implementation Status

### Current state of `dwba.cpp`
```
WavElj:          ✅ Fixed (inner zeroing, Wronskian matching, alpha normalization)
CalculateBoundState: ✅ Fixed (log-scale tracking for outward/inward matching)
InelDc:          ⚠️  Correct structure, but missing the FACTOR normalization
XSectn:          ⚠️  Missing FACTOR, wrong phase i^(La+Lo+2Lx+1), wrong 1/√(2La+1)
```

### What's left (the ~10^10 remaining error)
```
1. Missing FACTOR = 2*sqrt(ka*kb / (Ecm_a * Ecm_b)) ≈ 0.074
   → Apply in XSectn to the prefactor
   → Effect: reduces result by factor ~0.074² ≈ 0.005 (×200 reduction)

2. Wrong phase convention i^(La+Lo+2Lx+1) vs current implementation
   → Fixes sign/imaginary-unit mixing in coherent sum

3. Missing 1/√(2*La+1) normalization from SFROMI
   → La=2: factor of 1/√5 = 0.447 per amplitude

4. Possible additional factor from 4π / (ka*kb) in fLlsj definition (handbook eq.)
   vs our integral measure
```

---

## The Ten Commandments of Ptolemy

1. **Finite-range, not zero-range** — always 2D integral in (r_α, r_β)
2. **Prior form** — use V_bx (d-p potential) as the interaction W_β
3. **Stripping sign** — T2 = −S1 for the (d,p) mass factor
4. **Outward Numerov with inner zeroing** — for L>0 distorted waves
5. **Wronskian S-matrix extraction** — at r_max using Coulomb functions F_L, G_L
6. **alpha normalization** — chi = u * alpha at every grid point
7. **FACTOR = 2√(ki·ko/Ei·Eo)** — applied in SFROMI, not XSECTN
8. **Phase i^(La+Lo+2Lx+1)** — from partial wave expansion conventions
9. **1/√(2La+1)** normalization per distorted wave channel
10. **Units: ×10** for fm² → mb conversion

---

*Generated by Dudu 🐱 — 2026-03-12*  
*Source: `source.mor`, `source.f`, `handbook_of_direct_nuclear_reaction_for_retarded_theorist_v1-1.pdf`*
