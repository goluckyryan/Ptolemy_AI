# Ptolemy Inelastic Scattering: Theory & Implementation

**Dudu's guide for Ryan 🐱**
*Companion to `PTOLEMY_THEORY.md` (covers elastic + transfer DWBA)*
*Reaction: 206Hg(d,d')206Hg\*(4⁺, Ex=1.0 MeV, BELX=0.12, Elab=14.78 MeV)*
*Last updated: 2026-04-05*

---

## Table of Contents

1. [Overview](#1-overview)
2. [Physics: Inelastic vs Transfer](#2-physics-inelastic-vs-transfer)
3. [Input Parameters: BETA, BELX, and Deformation Lengths](#3-input-parameters-beta-belx-and-deformation-lengths)
4. [Form Factor: Nuclear and Coulomb](#4-form-factor-nuclear-and-coulomb)
5. [Inelastic Pipeline in Ptolemy](#5-inelastic-pipeline-in-ptolemy)
6. [INGRST / INRDIN: The Radial Integral](#6-ingrst--inrdin-the-radial-integral)
7. [BETCAL: Building the Amplitude](#7-betcal-building-the-amplitude)
8. [AMPCAL + XSECTN: DCS Assembly](#8-ampcal--xsectn-dcs-assembly)
9. [Quadrature and Grid Parameters](#9-quadrature-and-grid-parameters)
10. [WAVELJ for Inelastic](#10-wavelj-for-inelastic)
11. [C++ Implementation Notes](#11-c-implementation-notes)
12. [Known Issues and Bug Log](#12-known-issues-and-bug-log)
13. [Current Status](#13-current-status)

---

## 1. Overview

Ptolemy computes the **inelastic DWBA differential cross section** for collective excitation of the target.
The reaction is described as: projectile scatters off target, exciting it from ground state to state J\*.
No nucleon is transferred — the coupling is through the deformed nuclear potential.

Example: 206Hg(d,d')206Hg\*(4⁺, Ex=1.0 MeV)
- Projectile: deuteron (d), ejectile: deuteron (d)
- Target: 206Hg (ground state J=0), residual: 206Hg\* (J=4⁺, Ex=1.0 MeV)
- Coupling multipolarity: LX=4 (same as Jf for 0+ → 4+ transition)

The inelastic cross section is fundamentally different from transfer:
- **Transfer:** integrates over a 2D grid (r_α, r_β)
- **Inelastic:** integrates over a **1D grid** in r (the coupled radial integral via INRDIN/INGRST)

---

## 2. Physics: Inelastic vs Transfer

### 2.1 DWBA for Collective Excitation

In first-order DWBA, the inelastic S-matrix element (for entrance partial wave LI, exit LO) is:

$$S(L_I, L_O) = \frac{-i}{\hbar v} \int_0^\infty \chi^{(-)*}_{L_O}(r) \cdot V_{\text{coupl}}(r) \cdot \chi^{(+)}_{L_I}(r) \, dr$$

where:
- $\chi^{(+)}_{L_I}(r)$ = incoming distorted wave (entrance channel, partial wave LI)
- $\chi^{(-)*}_{L_O}(r)$ = outgoing distorted wave (exit channel, partial wave LO)
- $V_{\text{coupl}}(r)$ = coupling form factor (the derivative of the optical potential, weighted by deformation)

### 2.2 The Coupling Form Factor (Collective Model)

For a deformed nucleus with multipolarity LX, the nuclear coupling potential is:

$$V_{\text{coupl}}^{(N)}(r) = \beta_N \cdot R_N \cdot \frac{dU(r)}{dr}$$

where:
- $\beta_N$ = nuclear deformation parameter (dimensionless)
- $R_N = r_0 \cdot A^{1/3}$ = nuclear radius (fm)
- $\frac{dU}{dr}$ = radial derivative of the (complex) optical potential

For the Coulomb coupling (electric multipole LX):

$$V_{\text{coupl}}^{(C)}(r) = \begin{cases}
-\frac{3 Z_p Z_t e^2}{2L_X + 1} \cdot \frac{r^{L_X}}{R_C^{2L_X+1}} & r < R_C \\
-\frac{3 Z_p Z_t e^2}{2L_X + 1} \cdot \frac{R_C^{L_X}}{r^{L_X+1}} & r \geq R_C
\end{cases}$$

times $\beta_C$, with the coupling strength:

$$\text{VC} = -3 Z_p Z_t e^2$$

(Fortran source.f line 29461)

### 2.3 Total Coupling

The full form factor is:

$$F(r) = \beta_N R_N \frac{dU}{dr} + \beta_C \cdot V_C(r)$$

In practice, for 206Hg (heavy nucleus, Zp=1, Zt=80), the nuclear term dominates and the Coulomb coupling is a small correction.

---

## 3. Input Parameters: BETA, BELX, and Deformation Lengths

This is one of the most confusing parts of Ptolemy's input format.

### 3.1 Ptolemy Input Keywords

Ptolemy accepts **three ways** to specify coupling strength:

| Keyword | Meaning | Units |
|---------|---------|-------|
| `BETA`  | Nuclear deformation parameter β_N | dimensionless |
| `BETAC` | Coulomb deformation parameter β_C | dimensionless |
| `BELX`  | B(ELX) — reduced electric transition probability | e²·fm^(2LX) |

These are **not independent** — Ptolemy converts whichever you provide.

### 3.2 The BELX → β Conversion (MEBROT / BASCPL)

When `BELX` is given, Ptolemy computes β_C via the sharp-cutoff rotational model:

$$\beta_C = \frac{4\pi}{3 Z_t} \cdot \frac{\sqrt{B(EL_X)}}{R_C^{L_X}}$$

where $R_C = R_{C0} \cdot A^{1/3}$ is the Coulomb radius.

If `BETAC` is not given separately, Ptolemy sets:
$$\beta_N = \frac{R_C}{R_N} \cdot \beta_C \quad \text{(equal deformation lengths: } \delta_N = \delta_C\text{)}$$

### 3.3 Verified Values for 206Hg(d,d')206Hg*(4+)

With `BELX=0.12` (B(E4) = 0.12 e²·fm⁸), `RC0=1.303`, `A=206`, `J_f=4`:

```
R_C  = 1.303 × 206^(1/3) = 1.303 × 5.9025... = 7.692 fm
β_C  = 4π/(3×80) × sqrt(0.12) / 7.692^4 = 0.05172
R_N  = 1.151 × 206^(1/3) = 6.795 fm
β_N  = (R_C/R_N) × β_C = (7.692/6.795) × 0.05172 = 0.05855
```

> **Key fact:** `BELX` is B(E4), NOT β. If you put `BETA=0.12` you'd be setting β_N=0.12, which is ~2× too large and gives a cross section ~4× too big.

---

## 4. Form Factor: Nuclear and Coulomb

### 4.1 Nuclear Form Factor — Spline Derivative Approach

Ptolemy does **not** evaluate dU/dr analytically. Instead it:

1. Stores the optical potential V(r) on a uniform grid in Numerov-scaled form:
   ```
   POT[i] = -(H²/12E) × V(r_i)
   ```
   where `H = k × step` (dimensionless), `E = E_cm` (MeV), `step = 0.1` fm.

2. Fits a natural cubic spline (SPLNCB) to this scaled potential.

3. At each Gauss quadrature point r, evaluates the **spline derivative**:
   ```fortran
   XX = ALLOC(LB+N) + X*(2*ALLOC(LC+N) + X*3*ALLOC(LD+N))
   ```
   where LB/LC/LD are spline b/c/d coefficient arrays, X is fractional position in interval.

4. Recovers dV/dr via:
   ```fortran
   TERM = (12*E/H²) × R2S(1) × WT × XX
   ```
   where `R2S(1) = R_real` (fm), `WT` = quadrature weight.

The `(12E/H²)` factor converts the Numerov-scaled spline derivative back to physical units (MeV/fm).

> **Verified:** Spline derivative matches analytical dU/dr to 0.0001% (session 43).

### 4.2 Critical: Combined Imaginary Potential

The imaginary potential has **two terms** (volume + surface), but they must be combined into **one array** before spline:

```
V_imag_total(r) = V_vol × WS(r; ri, ai) + 4*V_surf × dWS/dr × ai
```

Then the spline is taken of `V_imag_total`, and the derivative multiplied by `R_imag` (NOT R_surf separately).

> **Bug fixed session 43:** Old code used separate splines for vol and surf imaginary terms, then multiplied each by their own radius. Since `R_vol ≠ R_surf`, this introduced up to 80% DCS error via interference. The fix: combine first, spline once, multiply by `R_imag` only.

### 4.3 BETRAT — The β Scaling Factor

The form factor is stored without β. The β is applied later in INRDIN:

```fortran
HR = HNUCR * ALLOC(LBNRAT + LX/2)
HI = HNUCI * ALLOC(LBNRAT + LX/2)
```

where `LBNRAT` stores `β_N × R_N` (the "deformation length" δ_N = β_N × R_N in fm).

In C++:
```cpp
double betaRat = betaN * R2S_real;   // δ_N = β_N × R_N (fm)
HR *= betaRat;
HI *= betaRat;
```

### 4.4 Coulomb Form Factor

For r ≥ R_C (exterior):
```
V_C(r) = VC / r^(LX+1)    where VC = -3 × Z_p × Z_t × e²  (in MeV·fm^LX)
```

For r < R_C (interior):
```
V_C(r) = VC × r^LX / R_C^(2LX+1)
```

In the Fortran INRDIN line ~22648:
```fortran
HCOUL = -WT * VC / r^(LX+1)    (exterior)
```

Note: `BETRAT` for Coulomb is `β_C × R_C` (or equivalently `β_C` when VC already contains the radius factors).

> **Bug fixed:** Early C++ had an extra `R_C^LX / (2LX+1)` factor → 390× too large for LX=4. Fixed by matching the Fortran VC definition exactly.


## 5. Inelastic Pipeline in Ptolemy

```
INPUT PARSING
    │  BELX → β_C → β_N (via MEBROT/BASCPL)
    ▼
ELASTIC WAVEFUNCTIONS (WAVELJ)
    │  Solve optical model Schrödinger for LI=0..Lmax (entrance channel)
    │  Solve optical model Schrödinger for LO=0..Lmax (exit channel)
    ▼
INGRST  ← Sets up 1D radial grid, potential arrays, quadrature
    │  Calls MAKPOT to build V(r) on grid
    │  Splines the potential → form factor storage
    ▼
INRDIN  ← Radial integral loop over (LI, LO, LX)
    │  For each (LI, LO): evaluates ∫ χ_LO*(r) × F(r) × χ_LI(r) dr
    │  Uses Gauss quadrature (CUBMAP) + 5-pt Lagrange interpolation of χ
    │  Multiplies by β × R (BETRAT)
    │  → S_matrix(LI, LO)
    ▼
BETCAL  ← Assembles partial-wave amplitudes
    │  For each (LO, MX): β(MX,LO) = (0.5/k_in) × Σ_LI (2LI+1)
    │    × CG(LI,0; LX,MX; LO,MX) × e^{i(σ_LI + σ_LO)} × S(LI,LO)/i
    ▼
AMPCAL  ← Angle-dependent amplitude
    │  F(MX, θ) = Σ_LO P_LO^MX(cos θ) × β(MX, LO)
    ▼
XSECTN  ← Differential cross section
    │  dσ/dΩ = 10 × Σ_MX w(MX) × |F(MX,θ)|²
    │  w(MX) = 1 for MX=0, 2 for MX>0
    │  Factor 10: fm² → mb/sr
    ▼
OUTPUT: dσ/dΩ (mb/sr) vs θ_CM
```


## 6. INGRST / INRDIN: The Radial Integral

### 6.1 INGRST — Grid Setup

INGRST (Fortran line ~21566) initializes the 1D radial grid for the inelastic integral:

- Grid: `r = step × i`, `i = 1..NSTEP`, `step = 0.1` fm
- `SUMMAX` = radial cutoff (fm) — default 20 fm for inelastic (not 35 like transfer)
- `NSTEP = SUMMAX / step`
- Grid extended by 20 extra steps past SUMMAX for spline edge stability

INGRST then calls:
1. `MAKPOT(NWP=3)` — builds the optical potential on the grid (both entrance and exit channels)
2. `SPLNCB` — fits cubic spline to the Numerov-scaled potential

### 6.2 INRDIN — The Integral Loop

INRDIN (Fortran line ~21894) evaluates:

$$S(L_I, L_O) = \sum_{\text{quad pts}} w_j \cdot \chi_{L_O}^*(r_j) \cdot F(r_j) \cdot \chi_{L_I}(r_j)$$

where the sum is over Gauss quadrature points in the interval [SUMMIN, SUMMAX].

**Key implementation details:**

1. **CUBMAP quadrature:** Maps [SUMMIN, SUMMAX] → [SUMMID, SUMMAX] with cubic density enhancement near SUMMID (where the integrand peaks). This is identical to the transfer CUBMAP but in 1D.

2. **Wavefunction at quadrature points:** WAVELJ is called with `NUMPTS = NUMPT` (~22 for 206Hg), which uses **5-point Lagrange interpolation** to evaluate χ(r) at the Gauss points. This is different from elastic (NUMPTS=0, uniform grid only).

3. **Form factor evaluation:** The spline derivative is evaluated at each Gauss point r_j using the stored spline coefficients.

4. **β application:** After the integral, multiply by `β_N × R_N` (nuclear) and `β_C` (Coulomb separately).

### 6.3 GAMSUM Default

For inelastic scattering, the Fortran default is:

```
GAMSUM = 5.0
```

(Not 1.0 as used in transfer.) GAMSUM controls the density of Gauss points near SUMMID. Getting this wrong causes systematic errors in the quadrature.

> **Bug fixed 2026-04-03:** GAMSUM was set to 1.0 (transfer default) instead of 5.0 (inelastic default). This was a ~5% error source.

### 6.4 Tail Correction (COULIN — in progress)

For high partial waves where the integrand extends beyond SUMMAX, Ptolemy applies an **asymptotic tail correction** using Coulomb functions:

$$I_{\text{tail}}(L_I, L_O) = \frac{1}{4} \cdot C \cdot \text{BETRAT} \cdot \left[(1+S_O)(1+S_I)FF - (1-S_O)(1-S_I)GG + i(1+S_O)(1-S_I)FG + i(1-S_O)(1+S_I)GF\right]$$

where FF/GG/FG/GF are integrals of Coulomb functions evaluated via RCASYM+CLINTS+COULIN.

**Status:** COULIN is ported but the tail correction currently makes DCS *worse* (19% mean vs 1.7% baseline) — indexing/sign bug under investigation.


## 7. BETCAL: Building the Amplitude

### 7.1 Formula

For each exit partial wave LO and magnetic quantum number MX:

$$\beta(M_X, L_O) = \frac{1}{2k_{\text{in}}} \sum_{L_I} (2L_I+1) \cdot C(L_I, 0; L_X, M_X | L_O, M_X) \cdot \frac{S(L_I, L_O)}{i} \cdot e^{i(\sigma_{L_I} + \sigma_{L_O})}$$

where:
- $k_{\text{in}}$ = entrance channel wave number (fm⁻¹)
- $C(\cdot|\cdot)$ = Clebsch-Gordan coefficient
- $S(L_I, L_O)$ = inelastic S-matrix element from INRDIN
- $\sigma_L = \arg\Gamma(L+1+i\eta)$ = Coulomb phase shift

### 7.2 Fortran Implementation (BETCAL lines 3620–3795)

```fortran
PHASE = SPHASE(KOFFS, I) + SIGIN(LI+1) + SIGOT(LO+1)
! SPHASE = nuclear phase of S-matrix element (arg of S/|S|)
! SIGIN(LI+1) = Coulomb phase σ_LI
! SIGOT(LO+1) = Coulomb phase σ_LO

SMATR = |S| × sin(PHASE)    ! Re(S/i) = Im(S)
SMATI = -|S| × cos(PHASE)   ! Im(S/i) = -Re(S)

! Accumulate into β(MX, LO) via Clebsch-Gordan × TEMPS
BETAS(MX, LO) += T(MX) × (SMATR + i×SMATI)
```

After accumulation, apply sqrt-factorial normalization for MX > 0:

```
β(MX, LO) *= ∏_{m=1}^{MX} 1/√((LO+m)(LO-m+1))
```

This normalizes the associated Legendre convention (no Condon-Shortley phase).

### 7.3 Seven Bugs Fixed in BETCAL

These were all verified during Sessions 17–25:

1. **CG loop:** BETCAL accumulates ALL (LI,LO) pairs into β(MX,LO) for each MX (not just diagonal)
2. **Associated Legendre:** Uses P_LO^MX **without** Condon-Shortley phase (matches Ptolemy's PLMSUB)
3. **Coulomb phases:** Must add σ_in(LI) + σ_out(LO) via SIGIN/SIGOT arrays
4. **Sqrt-factorial normalization:** Applied after accumulation, not inside loop
5. **FACTOR = 0.5/k_in:** (Not 1.0, not 0.5/k_out)
6. **CG argument:** CG(LI,0; LX,MX; LO,MX) — not CG(LX,MX; LI,0; LO,MX)
7. **MX weight:** w(MX)=2 for MX>0, w(0)=1 in DCS sum (factor FMNEG in Fortran)

> **Verified:** BETCAL produces 0.25% DCS error when fed the Fortran S-matrix directly (session 43). The remaining error is in the radial integral.


## 8. AMPCAL + XSECTN: DCS Assembly

### 8.1 AMPCAL

For each angle θ and each MX:

$$F(M_X, \theta) = \sum_{L_O} P_{L_O}^{M_X}(\cos\theta) \cdot \beta(M_X, L_O)$$

Uses associated Legendre polynomials P_L^M (Fortran: PLMSUB routine).
No Condon-Shortley phase (matches Ptolemy convention).

### 8.2 XSECTN — Differential Cross Section

$$\frac{d\sigma}{d\Omega}(\theta) = 10 \cdot \sum_{M_X=0}^{L_X} w(M_X) \cdot |F(M_X, \theta)|^2$$

where:
- Factor 10: converts fm² → mb/sr
- $w(0) = 1$, $w(M_X > 0) = 2$ (two orientations, M_X and -M_X)

### 8.3 Total Cross Section

Ptolemy also computes:

$$\sigma_{\text{total}} = \sum_{L_O} \frac{40\pi}{2L_O+1} \cdot \sum_{M_X} w(M_X) \cdot |\beta(M_X, L_O)|^2$$

(Fortran BETCAL: `SIGTOT += 40π/(2LO+1) × |β|²`, ×2 for MX>0)


## 9. Quadrature and Grid Parameters

### 9.1 Gauss Quadrature (CUBMAP)

The 1D integral over [SUMMIN, SUMMAX] uses Gauss-Legendre quadrature with a **cubic density mapping** (MAPSUM=2):

- Points concentrated near SUMMID (where the integrand peaks)
- GAMSUM controls the mapping concentration:
  - `GAMSUM = 5.0` (inelastic default, **not** 1.0)
  - Higher GAMSUM → more points near SUMMID

Number of quadrature points:

```
NUMPT = max(
  ceil((SUMMAX - SUMMIN) × SUMPTS × (k_in + k_out) / (4π)),
  NPSUM
)
```

where `SUMPTS=6`, `NPSUM=15`.

### 9.2 Verified Values for 206Hg (Elab=14.78 MeV)

| Parameter | Value |
|-----------|-------|
| SUMMIN | 0.0 fm |
| SUMMID | 10.0 fm |
| SUMMAX | 20.0 fm |
| NUMPT | 22 |
| MAPSUM | 2 (cubic) |
| GAMSUM | 5.0 ✅ (was 1.0, bug fixed 2026-04-03) |
| step | 0.1 fm |
| NSTEP | ~200 |

### 9.3 Wavefunction Interpolation

WAVELJ uses **5-point Lagrange interpolation** to evaluate χ(r) at Gauss points:

```
χ(r_gauss) ≈ Σ_{j=-2}^{+2} L_j(x) × χ(r_{i0+j})
```

where `i0` is the nearest grid index, `x` is the fractional offset.

This is controlled by `NUMPTS=NUMPT` in the WAVELJ call (vs `NUMPTS=0` for elastic/transfer which use only the uniform grid).


## 10. WAVELJ for Inelastic

### 10.1 Differences from Elastic / Transfer

| Mode | NUMPTS | Grid | Purpose |
|------|--------|------|---------|
| Elastic (WAVELJ) | 0 | Uniform only | S-matrix matching |
| Transfer (WAVELJ) | 0 | Uniform only | GRDSET uses own interpolation |
| **Inelastic (WAVELJ)** | **NUMPT (~22)** | Uniform + **5-pt Lagrange at Gauss pts** | Radial integral points |

For inelastic, WAVELJ stores the wavefunction at the uniform grid **and** evaluates it at `NUMPT` Gauss quadrature points via 5-point Lagrange, storing the results in a special output array.

### 10.2 Normalization Convention

WAVELJ normalizes so that at large r:

$$\chi_L(r) \xrightarrow{r\to\infty} \frac{i}{2}\left[H_L^-(kr) - S_L \cdot H_L^+(kr)\right]$$

where $H_L^\pm = G_L \pm iF_L$ are Coulomb-Hankel functions.

Matching code (Fortran lines 36029–36030):
```fortran
A1 = 0.5*(F*(1+SJR) + SJI*G)
A2 = 0.5*(G*(1-SJR) + SJI*F)
ALPHA = (WAV·A1 + WAV·A2) / |WAV|²
```

This convention is **identical** to the elastic solver — verified to match at r ≥ 18 fm (ratios = 1.0000).

### 10.3 S-Matrix Accuracy

| Component | Error vs Cleopatra 32-bit |
|-----------|--------------------------|
| Elastic S-matrix | 0.005% mean ✅ |
| Inelastic S-matrix (INRDIN) | 0.6% × (1.006 ± 0.030) — systematic bias |

The inelastic S-matrix has a mean ratio of 1.006 vs Fortran with σ=0.030. High-L tail elements (small magnitudes) have the largest fractional errors and dominate the 11.2% max DCS error.


## 11. C++ Implementation Notes

### 11.1 Key Files

| File | Purpose |
|------|---------|
| `test_inelastic.cpp` | Main driver: kinematics, WAVELJ calls, INRDIN, BETCAL, DCS |
| `src/dwba/ineldc.cpp` | INRDIN: radial integral loop |
| `src/dwba/grdset.cpp` | INGRST: grid setup, potential spline |
| `src/dwba/wavelj.cpp` | WAVELJ: Numerov solver + Lagrange interpolation |
| `src/dwba/xsectn.cpp` | BETCAL + AMPCAL + XSECTN: DCS assembly |
| `src/dwba/coulin.cpp` | RCASYM + CLINTS + COULIN: tail correction (in progress) |
| `include/coulin.h` | Headers for tail correction routines |

### 11.2 Build Command

```bash
cd ~/working/ptolemy_2019/Cpp_AI && g++ -std=c++17 -O2 -I include \
  test_inelastic.cpp \
  src/dwba/dwba.cpp src/dwba/bound.cpp src/dwba/coulin.cpp \
  src/dwba/rcwfn.cpp src/dwba/math_utils.cpp \
  src/dwba/wavelj.cpp src/dwba/setup.cpp src/dwba/potential_eval.cpp \
  src/dwba/spline.cpp src/dwba/xsectn.cpp \
  src/dwba/grdset.cpp src/dwba/ineldc.cpp \
  src/dwba/a12.cpp src/dwba/av18_potential.cpp \
  src/input/Isotope.cpp \
  /tmp/stub_zr.cpp \
  -o test_inelastic_coulin -lm
```

### 11.3 Fortran Reference Input

```
REACTION: 206Hg(d,d)206Hg(4+1.000) ELAB=14.780
 INCOMING
V=96.891 r0=1.151 a=0.793 VI=2.023 ri0=1.322 ai=0.264 VSI=10.378 rsi0=1.360 asi=0.897 rc0=1.303
;
 OUTGOING
V=96.891 r0=1.151 a=0.793 VI=2.023 ri0=1.322 ai=0.264 VSI=10.378 rsi0=1.360 asi=0.897 rc0=1.303
;
JBIGA=0
BETA=0.12
PRINT 2
ANGLEMIN=0 ANGLEMAX=180 ANGLESTEP=1
;
```

> Note: use `(d,d)` not `(d,d')` — Ptolemy cannot parse the prime character.

### 11.4 Kinematics (206Hg, Elab=14.78 MeV)

```
m_d  = 2.01410 u = 1875.61 MeV/c²
M_Hg = 205.97 u  = 191840.7 MeV/c²
μ_in = m_d × M_Hg / (m_d + M_Hg) = 1869.8 MeV/c²
E_cm = Elab × M_Hg/(m_d + M_Hg) = 14.78 × 191840.7/193716.3 = 14.637 MeV
k_in = √(2μ E_cm) / ħc = √(2 × 1869.8 × 14.637) / 197.327 = 1.4779 fm⁻¹
η_in = Z_p × Z_t × e² × μ / (ħ² k_in) = 1×80×1.44 MeV·fm × 1869.8/(197.327² × 1.4779) = 27.9

(exit channel same as entrance for inelastic — same projectile, same target up to Q)
Q = -Ex = -1.0 MeV (excited state at 1.0 MeV)
E_cm_out = E_cm_in + Q = 14.637 - 1.0 = 13.637 MeV
k_out = 1.4262 fm⁻¹
```


## 12. Known Issues and Bug Log

### 12.1 Bugs Fixed (chronological)

| Session | Bug | Impact |
|---------|-----|--------|
| 17–25 | 7 BETCAL bugs (CG, phases, normalization, etc.) | Large — DCS wrong |
| 42 | Coulomb FF: extra `R_C^LX/(2LX+1)` factor | ~390× for LX=4 (but nuclear dominates) |
| 43 | R_imag: separate vol/surf splines → combine first | Up to 80% via interference |
| 43 | INGRST grid: extend 20 steps past SUMMAX for spline | Edge instability |
| 46 | ElasticSolver: one WS for both real+imag → separate calls | DCS: 15% → 10% |
| 47 | Replaced ElasticSolver with WAVELJ directly | 10% → 1.7% mean |
| 48 | BELX vs BETA confusion (BELX=B(E4), not β_N) | ~4× DCS |
| 2026-04-03 | GAMSUM=1.0 → 5.0 (inelastic default) | ~5% |
| 2026-04-03 | Wavefunction interpolation: indices wf[i0-2..i0+2] not wf[i0-1..i0+3] | ~2% |

### 12.2 Active Issue: COULIN Tail Correction

**Status:** Ported but incorrect — makes DCS worse.

Without correction: 1.72% mean, 11.2% max
With correction: 19.48% mean, 54.6% max

Likely culprits:
1. COULIN array indexing (Fortran 1-based → C++ 0-based)
2. Sign or phase in IRTOIN formula
3. IFACTR (imaginary unit i factor) treatment
4. Whether Fortran reference data already includes the tail correction

### 12.3 Remaining DCS Error

| Metric | Value |
|--------|-------|
| Mean DCS error | 1.72% |
| Max DCS error | 11.2% (at θ~24°, from high-L tail S-matrix elements) |
| S-matrix ratio (mean) | 1.006 ± 0.030 vs Fortran |
| BETCAL check (Fortran S input) | 0.25% ✅ |
| Elastic S-matrix | 0.005% ✅ |

The error is localized to high-L S-matrix elements with small magnitudes — likely the COULIN tail correction or the Lagrange interpolation at extreme radii.

### 12.4 What Has Been Ruled Out

- BETCAL formula ✅ (verified with Fortran S-matrix)
- Coulomb form factor ✅ (fixed, minor impact)
- Spline vs analytical dV/dr ✅ (0.0001% difference)
- CUBMAP Gauss quadrature mapping ✅
- Coulomb phases ✅
- Spin CG coefficients ✅ (= 1.0 trivially, J_in=0, J_out=8)
- Kinematics (k, η, μ) ✅ (verified vs Fortran)
- WAVELJ wavefunction normalization ✅


## 13. Current Status

**As of 2026-04-05:**

### ✅ Complete
- Elastic S-matrix: 0.005% match vs Cleopatra 32-bit
- BETCAL: 0.25% DCS error with Fortran S-matrix input
- Form factor (dV/dr via spline): 0.0001% vs analytical
- Kinematics: verified
- WAVELJ integration: working (uses inelastic NUMPTS mode)

### 🔄 In Progress
- COULIN tail correction: ported (RCASYM ✅, CLINTS ✅, COULIN ✅) but result is wrong
  - Without: 1.72% mean DCS error
  - With: 19.48% mean (worse — indexing/sign bug)

### 📋 Next Steps
1. Debug COULIN: print actual FF/FG/GF/GG values vs Fortran PRINT=5 output
2. Verify ID/IL index mapping for specific (LI,LO) pairs
3. Check whether Fortran reference DCS already includes tail correction
4. Once COULIN fixed: re-validate full pipeline

### 🎯 Target
- Match Cleopatra 32-bit to <0.1% mean DCS (like transfer reaction)
- Current: 1.72% mean without tail correction

---

*Reference binary:* `~/working/ptolemy_2019/digios/analysis/Cleopatra/ptolemy` (32-bit, INELOCA1, no SO)
*Reference data:* `Cpp_AI/fortran_inel_ref.dat`
*Standing rules:* See `ptolemy_rules.md` — 32-bit only, 206Hg only, no Fortran compilation

