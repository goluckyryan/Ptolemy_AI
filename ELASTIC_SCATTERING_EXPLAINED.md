# How Ptolemy Computes Elastic Scattering

**Dudu's reference for Ryan 🐱**  
*Cross-referenced with: `handbook_of_direct_nuclear_reaction_for_retarded_theorist_v1-1.pdf` (§2, §4) and `source.mor`*

---

## The Full Flow (Elastic Only)

```
INPUT: optical potential parameters + kinematics
    │
    ▼
WAVPOT  ← build potential arrays V(r) on radial grid
    │
    ▼
WAVEF / WAVELJ  ← solve optical model ODE for each L, extract S_LJ
    │
    ▼
BETCAL  ← compute angle-independent β coefficients from S-matrix
    │
    ▼
AMPCAL  ← for each angle θ, sum β·P_L(cosθ) → amplitude f(θ)
    │        also compute Coulomb amplitude f_C(θ) analytically
    │
    ▼
ELDCS   ← |f(θ) + f_C(θ)|² × 10 → dσ/dΩ in mb/sr
```

---

## Step 1: The Schrödinger Equation Being Solved

For each partial wave L (and total J = L ± ½ if spin-orbit is on), Ptolemy solves:

$$\left[\frac{d^2}{dr^2} + k^2 - \frac{L(L+1)}{r^2} - \frac{2\mu}{\hbar^2}\bigl(V_C(r) + U_{opt}(r)\bigr)\right]\chi_{LJ}(r) = 0$$

where:
- $k^2 = 2\mu E_{cm}/\hbar^2$ (wave number squared, real)
- $V_C(r)$ = Coulomb potential (point charge, real)
- $U_{opt}(r)$ = complex optical potential (nuclear — absorptive imaginary part)
- $\Lambda(L,s,J) = \frac{1}{2}[J(J+1) - L(L+1) - s(s+1)]$ is the spin-orbit coupling factor

**Normalization (boundary condition at large r):**

$$\chi_{LJ}(r) \xrightarrow{r\to\infty} \frac{e^{i\sigma_L}}{2i}\left[S_{LJ}H_L^+(kr) - H_L^-(kr)\right]$$

where $\sigma_L = \arg\Gamma(L+1+i\eta)$ is the Coulomb phase shift, $H_L^\pm = G_L \pm iF_L$, and $|S_{LJ}| \leq 1$ (= 1 for pure Coulomb).

**This is your Handbook eq. (27) exactly.**

---

## Step 2: WAVELJ — Numerical Solution (Outward Numerov)

```
Subroutine WAVEF (entry WAVELJ) in source.mor ~line 29969
```

### Numerov recursion
$$u_{i+1} = \frac{2\left(1 - \frac{5h^2}{12}f_i\right)u_i - \left(1 + \frac{h^2}{12}f_{i-1}\right)u_{i-1}}{1 + \frac{h^2}{12}f_{i+1}}$$

with $f_i = k^2 - L(L+1)/r_i^2 - (2\mu/\hbar^2)U(r_i)$.

Initial conditions: $u_0 = 0$, $u_1 = h^{L+1}$ (regular at origin, dominant power-law behavior $\sim r^{L+1}$).

### Overflow control for large L (the hard part)
For L ≥ 2, the centrifugal term $L(L+1)/r^2$ is huge at small $r$, making $f_i \ll 0$ and causing exponential growth:

- At r = 0.1 fm, L=2: $f \approx -600\,\text{fm}^{-2}$, growth per step $\approx e^{\sqrt{600}\times 0.1} \approx 11\times$
- After ~36 steps the solution exceeds `BIGNUM = 1e30`

**Ptolemy's fix** (source.mor ~30470): when $|u_{i+1}| >$ `BIGNUM`:
1. Divide all stored values by $|u_{i+1}|$ (rescale)
2. **Zero out the inner region** where $|u_j| <$ `STEPI` × $|u_{i+1}|$ after rescaling
3. Update `ISTRT` (first nonzero index) accordingly

This works because the physical regular solution IS zero in the deep interior (it only looks large due to numerical overflow of the power-law growth; after normalization via $\alpha$, the interior is unphysical).

### S-matrix extraction (Wronskian matching, source.mor ~30870)
At the outer boundary $r_N$ and $r_{N-1}$, match the Numerov solution to the asymptotic Coulomb functions. Solving the 2×2 linear system:

$$u_N = \alpha_R F_L(kr_N) + \alpha_I G_L(kr_N) + S_{LJ}\bigl(\alpha_I F_L - \alpha_R G_L\bigr)_N$$
$$u_{N-1} = \alpha_R F_L(kr_{N-1}) + \alpha_I G_L(kr_{N-1}) + S_{LJ}\bigl(\alpha_I F_L - \alpha_R G_L\bigr)_{N-1}$$

Subtracting and forming the Wronskian gives $S_{LJ}$ directly (no approximation).

### Normalization (source.mor ~30913–30916)
Once $S_{LJ}$ is known, back-calculate:

```fortran
A1 = 0.5*(F*(1+SJR) + SJI*G)    ! Re part of asymptotic chi / u_N
A2 = 0.5*(G*(1-SJR) + SJI*F)    ! Im part
ALPHAR = (WAVR(N)*A1 + WAVI(N)*A2) / (WAVR(N)^2 + WAVI(N)^2)
ALPHAI = (WAVR(N)*A2 - WAVI(N)*A1) / (WAVR(N)^2 + WAVI(N)^2)
```

Then `chi(r_i) = (ALPHAR + i*ALPHAI) * u_i` for all stored points.

**After this step**, `chi_LJ(r)` is the correctly normalized distorted wave, satisfying:
- `chi(0) = 0`
- `chi(r → ∞) → (e^iσ_L / 2i)(S_LJ H⁺ - H⁻)`

---

## Step 3: BETCAL — Angle-Independent β Coefficients

```
Subroutine BETCAL in source.mor line 3358
```

For elastic scattering (`ELSW = .TRUE.`), for each $(L_O, L_I, L_X)$ combination:

### S-matrix with Coulomb phases
```fortran
PHASE = SIGIN(LI+1) + SIGIN(LO+1)    ! σ_Li + σ_Lo
TR = cos(PHASE)
TI = sin(PHASE)
SR = S(1,KOFFS,I)    ! Re(S_LJ)
SI = S(2,KOFFS,I)    ! Im(S_LJ)
SMATR = SR*TI + SI*TR    ! Re(e^{2iσ} S - 1) piece
SMATI = SI*TI - SR*TR    ! Im(e^{2iσ} S - 1) piece
```

For $L_X = 0$ (elastic, no spin-orbit): this computes $e^{2i\sigma_L}(S_L - 1)$.

### Clebsch-Gordan coupling into TEMPS
```fortran
TEMPS(I+MX) = FACTOR * (2*LI+1) * CG(2*LI, 2*LX, 0, 2*MX, 2*LO, 2*MX)
```
where `FACTOR = 0.5/UK`. This is the coupling factor for each $(L_O, L_X, M_X, L_I)$.

### Accumulate β
```fortran
BETAS(1,KOFFZ+MX,ILO) += TEMPS(I+MX) * SMATR
BETAS(2,KOFFZ+MX,ILO) += TEMPS(I+MX) * SMATI
```

**What β represents:**
$$\beta_{L_X, M_X}(L_O) = \frac{1}{2k_\alpha}\sum_{L_I}(2L_I+1)\,C(L_I,0,L_X,M_X;L_O,M_X)\,e^{2i\sigma_{L_I}}(S_{L_I}-\delta_{L_X,0})$$

This is the coefficient of $P_{L_O}^{M_X}(\cos\theta)$ in the nuclear scattering amplitude.

---

## Step 4: AMPCAL — The Full Amplitude f(θ)

```
Subroutine AMPCAL in source.mor line 220
```

### Coulomb amplitude (analytic, source.mor ~345)
$$f_C(\theta) = -\frac{\eta}{2k\sin^2(\theta/2)}\exp\left[2i\sigma_0 - 2i\eta\ln\sin(\theta/2)\right]$$

This is the standard Rutherford scattering amplitude (your Handbook §2.3, Eq. 10).

```fortran
PHASE = 2*(SIGZRO - ETA*LOG(SIN(θ/2)))
TEMP  = -0.5*ETA / (UK * SIN²(θ/2))
FTC(1) = TEMP * cos(PHASE)    ! Re f_C
FTC(2) = TEMP * sin(PHASE)    ! Im f_C
```

### Nuclear amplitude (partial wave sum)
For each angular momentum $(L_O, L_X, M_X)$:
```
PLM_factor = P_{LO}^{MX}(cosθ)     ← associated Legendre function
```
Accumulate:
```fortran
F(1,KOFFS) += PLM * BETAS(1,KOFFS, LO-LBASE+1)    ! Re f_nuclear
F(2,KOFFS) += PLM * BETAS(2,KOFFS, LO-LBASE+1)    ! Im f_nuclear
```

### Final amplitude (elastic)
$$f(\theta) = f_C(\theta) + f_{\text{nuclear}}(\theta)$$

```fortran
F(IRI,1) = F(IRI,1) + FTC(IRI)   ! add Coulomb to nuclear amplitude
```

---

## Step 5: ELDCS — Cross Section

```
Subroutine ELDCS in source.mor line 12989
```

$$\frac{d\sigma}{d\Omega} = |f(\theta)|^2 \times 10 \quad \text{[mb/sr]}$$

(The ×10 converts fm² → mb.)

```fortran
SIGMA = 10 * (F(1,1)**2 + F(2,1)**2)    ! no spin-orbit
RUTH  = 10 * (FCOUL(1,3)**2 + FCOUL(2,3)**2)
SIGMAR = SIGMA / RUTH    ! ratio to Rutherford
```

---

## Summary: The Complete Elastic Formula

Putting it all together, Ptolemy computes (your Handbook §2.4, combining Eqs. 7 and 10):

$$\boxed{\frac{d\sigma}{d\Omega}(\theta) = \left|f_C(\theta) + f_N(\theta)\right|^2}$$

where:

$$f_C(\theta) = -\frac{\eta}{2k\sin^2(\theta/2)} e^{i[2\sigma_0 - 2\eta\ln\sin(\theta/2)]}$$

$$f_N(\theta) = \frac{1}{2ik}\sum_{L=0}^{L_{\max}}(2L+1)\,e^{2i\sigma_L}(S_L - 1)\,P_L(\cos\theta)$$

(For spin-orbit: $S_L \to S_{LJ}$, and there are separate amplitudes for each $J = L \pm \frac{1}{2}$.)

---

## Comparison: Your Handbook vs. Ptolemy

| Step | Your Handbook | Ptolemy |
|---|---|---|
| **ODE** | Eq. 5: $(-d²/dr² + l(l+1)/r²)u + 2μV u/ℏ² = k²u$ | Same in WAVELJ |
| **Regular BC at origin** | $u \sim r^{l+1}$ | `u[1] = h^{L+1}`, Numerov forward |
| **Asymptotic BC** | Eq. 27: $\chi → \frac{e^{iσ}}{2i}(SH^+ - H^-)$ | Same; ALPHAR/ALPHAI enforce this |
| **S-matrix extraction** | §2.7: integration formula (Wronskian) | Two-point Wronskian at $r_N$, $r_{N-1}$ |
| **Coulomb phase σ_L** | §2.3: $\arg\Gamma(L+1+iη)$ | Precomputed in `/cnstnt/SIGZRO` array |
| **Coulomb amplitude** | §2.3, Eq. 10 | `AMPCAL` lines ~345–390, exact formula |
| **Nuclear amplitude** | §2.4: $f_N = \frac{1}{2ik}\sum(2L+1)e^{2iσ}(S-1)P_L$ | `BETCAL` + `AMPCAL`: β×P_L sum |
| **Cross section** | §2.5: $\|f_C + f_N\|^2$ | `ELDCS`: `(F_r² + F_i²)*10` |
| **Spin-orbit** | §2.9: separate $J = L ± \frac{1}{2}$ amplitudes | `NSPL=2` path in `BETCAL`/`AMPCAL` |
| **Units** | fm²/sr | ×10 for mb/sr |

---

## The One Key Formula to Check in Our C++ Code

The elastic scattering amplitude is **exactly:**

$$f_N(\theta) = \frac{1}{2ik_\alpha}\sum_{L=0}^{L_{\max}}(2L+1)\,e^{2i\sigma_L}\,(S_L - 1)\,P_L(\cos\theta)$$

This requires:
1. ✅ **Correct $S_L$** — extracted by Wronskian from WAVELJ
2. ✅ **Correct Coulomb phases $\sigma_L$** — from `RCWFN`/`Rcwfn` for each L
3. ✅ **Correct $P_L(\cos\theta)$** — from `std::legendre` (C++17)
4. ✅ **Factor $(2L+1)/(2ik)$** — the standard partial wave prefactor
5. ✅ **Add $f_C$** — Rutherford amplitude on top of nuclear amplitude

---

## What to Validate First

Before touching the transfer reaction code, verify elastic scattering independently:

### Test: 33Si + d at 20 MeV, elastic channel
Expected: $S_L \approx e^{2i\delta_L}$ with $|S_L| < 1$ due to absorption.

**Quick checks:**
```
1. Plot |S_L| vs L — should be ~0 for L < L_grazing, ~1 for L > L_grazing
2. At θ = 0°: f_N(0) = (1/2ik) Σ(2L+1)(S_L-1) — should be mostly imaginary
3. Total reaction cross section: σ_R = (π/k²) Σ(2L+1)(1-|S_L|²) — should be ~few hundred mb
4. dσ/dΩ at small θ should → Rutherford
5. Compare our S_L with Ptolemy output (S-matrix is printed in Ptolemy output with IPRINT≥2)
```

### Reading Ptolemy's S-matrix output
In the Ptolemy output file, look for lines like:
```
L   REAL       IMAG       MAGNITUDE   PHASE
 0  -0.2962  -1.0213     1.064       ...
 1  ...
```
The `S_0 = (-0.2962, -1.0213)` that appears in our diagnostic IS the elastic S-matrix for L=0 — that's our normalization reference point.

---

## Concrete Diagnosis Plan

**Step 1:** Add a diagnostic that prints $S_L$ for all L (both channels: d+33Si and p+34Si).  
**Step 2:** Compare with Ptolemy's printed S-matrix (run with IPRINT=2 or higher).  
**Step 3:** Compute elastic cross section from our $S_L$ and verify against Ptolemy.  
**Step 4:** Only then move to the transfer reaction.

The elastic cross section is:
$$\frac{d\sigma}{d\Omega} = \left|f_C(\theta) + \frac{1}{2ik}\sum_{L}(2L+1)e^{2i\sigma_L}(S_L-1)P_L(\cos\theta)\right|^2 \times 10 \quad \text{mb/sr}$$

If our $S_L$ matches Ptolemy → our distorted waves are correct → the transfer issue is purely in the prefactor/coupling.

---

*Generated by Dudu 🐱 — 2026-03-12*  
*Sources: `source.mor` (BETCAL line 3358, AMPCAL line 220, ELDCS line 12989, WAVELJ ~30793), handbook §2, §4*
