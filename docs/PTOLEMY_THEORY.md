# Ptolemy DWBA: Theory & Algorithms

**Dudu's guide for Ryan 🐱**
*Cross-referenced with: `handbook_of_direct_nuclear_reaction_for_retarded_theorist_v1-1.pdf` and `source.mor`*
*Last updated: 2026-03-29*

---

## Table of Contents

1. [Overview](#1-overview)
2. [Elastic Scattering](#2-elastic-scattering)
3. [Bound States](#3-bound-states)
4. [DWBA Transfer Reaction](#4-dwba-transfer-reaction)
5. [Zero-Range vs Finite-Range](#5-zero-range-vs-finite-range)
6. [Coordinate Geometry](#6-coordinate-geometry)
7. [Integration Grid (GRDSET)](#7-integration-grid-grdset)
8. [Radial Integrals (INELDC)](#8-radial-integrals-ineldc)
9. [Transfer S-Matrix (SFROMI)](#9-transfer-s-matrix-sfromi)
10. [Cross Section (BETCAL → AMPCAL → XSECTN)](#10-cross-section)
11. [Comparison: Handbook vs Ptolemy](#11-comparison-handbook-vs-ptolemy)

---

## 1. Overview

Ptolemy computes the **finite-range DWBA differential cross section** for a transfer reaction A(a,b)B.

Example: 16O(d,p)17O — projectile a=d, ejectile b=p, transferred particle x=n, target A=16O, residual B=17O.

```
INPUT PARSING
    │
    ▼
BOUND STATES  ← solve Schrödinger eq. for n inside B and inside a
    │
    ▼
DISTORTED WAVES  ← solve optical model Schrödinger eq. for a+A and b+B
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

## 2. Elastic Scattering

Before calculating transfer reactions, we must describe the elastic scattering of the projectile and ejectile in the entrance and exit channels. This provides the "distorted waves" (χ).

### 2.1 The Schrödinger Equation

For each partial wave L (and total J = L ± ½ if spin-orbit is on), Ptolemy solves:

$$\left[\frac{d^2}{dr^2} + k^2 - \frac{L(L+1)}{r^2} - \frac{2\mu}{\hbar^2}\bigl(V_C(r) + U_{opt}(r)\bigr)\right]\chi_{LJ}(r) = 0$$

where:
- $k^2 = 2\mu E_{cm}/\hbar^2$ (wave number squared)
- $V_C(r)$ = Coulomb potential (point charge)
- $U_{opt}(r)$ = complex optical potential (nuclear + absorptive imaginary)
- $\Lambda(L,s,J) = \frac{1}{2}[J(J+1) - L(L+1) - s(s+1)]$ is the spin-orbit coupling factor

**Normalization (boundary condition at large r):**

$$\chi_{LJ}(r) \xrightarrow{r\to\infty} \frac{e^{i\sigma_L}}{2i}\left[S_{LJ}H_L^+(kr) - H_L^-(kr)\right]$$

where $\sigma_L = \arg\Gamma(L+1+i\eta)$ is the Coulomb phase shift, $H_L^\pm = G_L \pm iF_L$, and $|S_{LJ}| \leq 1$ (= 1 for pure Coulomb).

**(Handbook eq. 27)**

### 2.2 WAVELJ — Numerov Integration

```
Subroutine WAVEF (entry WAVELJ) in source.mor ~line 29969
```

**Numerov recursion:**
$$u_{i+1} = \frac{2\left(1 - \frac{5h^2}{12}f_i\right)u_i - \left(1 + \frac{h^2}{12}f_{i-1}\right)u_{i-1}}{1 + \frac{h^2}{12}f_{i+1}}$$

with $f_i = k^2 - L(L+1)/r_i^2 - (2\mu/\hbar^2)U(r_i)$.

Initial conditions: $u_0 = 0$, $u_1 = h^{L+1}$ (regular at origin).

#### Overflow control for large L

For L ≥ 2, the centrifugal term $L(L+1)/r^2$ causes exponential growth at small r:

- At r = 0.1 fm, L=2: $f \approx -600\,\text{fm}^{-2}$, growth per step $\approx e^{\sqrt{600}\times 0.1} \approx 11\times$
- After ~36 steps the solution exceeds `BIGNUM = 1e30`

**Ptolemy's fix** (source.mor ~30470): when $|u_{i+1}| >$ `BIGNUM`:
1. Divide all stored values by $|u_{i+1}|$ (rescale)
2. Zero out the inner region where $|u_j| <$ `STEPI` × $|u_{i+1}|$ after rescaling
3. Update `ISTRT` (first nonzero index)

This works because the physical regular solution IS zero in the deep interior — the large values are numerical artifacts of the power-law growth.

### 2.3 S-Matrix Extraction (Wronskian)

At the outer boundary $r_N$ and $r_{N-1}$, match the Numerov solution to asymptotic Coulomb functions:

$$S_L = \frac{u_{N-1} \cdot H^{+}_{L}(r_N) - u_N \cdot H^{+}_{L}(r_{N-1})}{u_{N-1} \cdot H^{-}_{L}(r_N) - u_N \cdot H^{-}_{L}(r_{N-1})}$$

where $H^\pm_L = G_L \pm iF_L$ and $F_L, G_L$ from RCWFN (Steed's method). This is exact — no approximation.

### 2.4 Normalization

Once $S_{LJ}$ is known, back-calculate:

```fortran
A1 = 0.5*(F*(1+SJR) + SJI*G)
A2 = 0.5*(G*(1-SJR) + SJI*F)
ALPHAR = (WAVR(N)*A1 + WAVI(N)*A2) / (WAVR(N)^2 + WAVI(N)^2)
ALPHAI = (WAVR(N)*A2 - WAVI(N)*A1) / (WAVR(N)^2 + WAVI(N)^2)
```

Then `chi(r_i) = (ALPHAR + i*ALPHAI) * u_i` for all stored points.

### 2.5 Elastic Cross Section

The complete elastic formula (Handbook §2.4):

$$\boxed{\frac{d\sigma}{d\Omega}(\theta) = \left|f_C(\theta) + f_N(\theta)\right|^2}$$

where:

$$f_C(\theta) = -\frac{\eta}{2k\sin^2(\theta/2)} e^{i[2\sigma_0 - 2\eta\ln\sin(\theta/2)]}$$

$$f_N(\theta) = \frac{1}{2ik}\sum_{L=0}^{L_{\max}}(2L+1)\,e^{2i\sigma_L}(S_L - 1)\,P_L(\cos\theta)$$

(For spin-orbit: $S_L \to S_{LJ}$, separate amplitudes for each $J = L \pm \frac{1}{2}$.)

### 2.6 Elastic Flow in Ptolemy

```
WAVPOT  → build potential arrays V(r) on radial grid
WAVELJ  → solve ODE for each L, extract S_LJ
BETCAL  → β coefficients from S-matrix: β(Lo) = (1/2k) Σ (2Li+1) e^{2iσ} (S-1) × CG
AMPCAL  → f(θ) = Σ β·P_L(cosθ) + f_C(θ)
ELDCS   → |f(θ)|² × 10 → dσ/dΩ in mb/sr
```

---

## 3. Bound States

### 3.1 The Problem

Solve the radial Schrödinger equation for a bound state:

$$\left[\frac{d^2}{dr^2} - \frac{l(l+1)}{r^2} + \frac{2\mu}{\hbar^2}\left(E - V_{WS}(r) - V_{SO}(r)\right)\right] u_{lj}(r) = 0$$

with $E < 0$ (binding energy), boundary condition $u(0) = 0$, and exponential decay at large r.

**Two bound states are needed for transfer:**
- $\phi_T(r)$ = transferred particle bound in residual B (e.g., n in 17O)
- $\phi_P(r)$ = transferred particle bound in projectile a (e.g., n in deuteron)

**Normalization:** $\int_0^\infty u_{lj}^2(r)\, dr = 1$, then $\phi_{lj}(r) = u_{lj}(r)/r$

### 3.2 Algorithm: Inward-Outward Matching

Unlike simple shooting methods (unstable for bound states), Ptolemy integrates from both sides and matches at $R_{match}$.

**Outward integration** ($0 \to R_{match}$):
- BC: $u(0) = 0$, $u(h) \sim h^{l+1}$
- Method: Numerov integration
- Counts nodes during integration

**Inward integration** ($R_{max} \to R_{match}$):
- BC: $u(R_{max}) \sim \exp(-\kappa R_{max})$ where $\kappa = \sqrt{2\mu |E|} / \hbar$
- Method: Numerov integration (backward)

**Matching condition** at $R_{match}$:
$$\text{Diff} = \left( \frac{u'}{u} \right)_{out} - \left( \frac{u'}{u} \right)_{in} \to 0$$

### 3.3 Depth Search Iteration

Iterate on potential depth $V$ to drive $\text{Diff} \to 0$:

1. **Node check**: If nodes > n, V too deep (decrease); if nodes < n, V too shallow (increase)
2. **Newton-Raphson**: Once node count is correct, update $V_{new} = V_{old} + \text{factor} \times \text{Diff}$
3. **Adaptive factor**: Starts at 2.0; halved when oscillation detected (sign of Diff flips)

---

## 4. DWBA Transfer Reaction

### 4.1 The Transition Amplitude

For a transfer reaction $a (= b + x) + A \to b + B (= A + x)$:

$$T_{ba} = \int d^3 r_b \int d^3 r_a \, \chi_{b}^{(-)\ast}(\vec{r}_{b}) \langle \psi_B \psi_b | V | \psi_A \psi_a \rangle \chi_{a}^{(+)}(\vec{r}_{a})$$

Separating internal states from radial motion:

$$T_{ba} \propto \int d^3 r_b \int d^3 r_a \, \chi_{b}^{\ast}(\vec{r}_{b}) \, \phi_{Bx}^{\ast}(\vec{r}_{x}) \, V(\vec{r}_{bx}) \, \phi_{ax}(\vec{r}_{bx}) \, \chi_{a}(\vec{r}_{a})$$

- $\phi_{Bx}$: Bound state wavefunction of x in B (target bound state)
- $\phi_{ax}$: Bound state wavefunction of x in a (projectile bound state)

### 4.2 Partial Wave Expansion

Expanding all wavefunctions in partial waves:

$$T \sim \sum_{L_a, L_b} \int r_{a}^{2} \, dr_a \int r_{b}^{2} \, dr_b \, \chi_{L_b}(r_b) \, K_{L_a L_b}(r_a, r_b) \, \chi_{L_a}(r_a)$$

The kernel $K$ involves angular integration:

$$K(r_a, r_b) = \int_{-1}^{1} d(\cos \theta) \, \phi_{Bx}^{\ast}(r_x) \, V(r_{bx}) \, \phi_{ax}(r_{bx}) \, P_L(\cos \theta)$$

where $r_x$ and $r_{bx}$ depend on $r_a, r_b, \theta$ via the law of cosines.

---

## 5. Zero-Range vs Finite-Range

### Zero-Range Approximation

Assumes the interaction and projectile bound state are pointlike:

$$V(\vec{r}_{bx}) \, \phi_{ax}(\vec{r}_{bx}) \approx D_0 \, \delta(\vec{r}_{bx})$$

This collapses the 6D integral to a 3D integral:

$$T_{\text{ZR}} \propto \int d^3 r \, \chi_{b}^{\ast}(\beta \vec{r}) \, \phi_{Bx}^{\ast}(\vec{r}) \, \chi_{a}(\vec{r})$$

Fast but inaccurate for heavy ion reactions, high energies, and states where finite size matters.

### Finite-Range (Ptolemy)

Ptolemy performs the full **finite-range** calculation — double radial integral + angular integral. Key differences:

| Aspect | Zero-Range | Finite-Range (Ptolemy) |
|---|---|---|
| Integral | 1D radial | 2D radial + angular |
| Projectile structure | Point (D₀δ) | Full φ_P(r) wavefunction |
| Computational cost | $N_r$ | $N_r^2 \times N_\theta$ |
| Recoil effects | Approximate | Exact |
| Accuracy | 10-30% off for (d,p) | Correct |

---

## 6. Coordinate Geometry

### 6.1 Jacobi Coordinates

For a transfer reaction, the relevant coordinates are:
- $\vec{r}_\alpha$ = relative coordinate (entrance channel, e.g., d ↔ 16O)
- $\vec{r}_\beta$ = relative coordinate (exit channel, e.g., p ↔ 17O)
- $\vec{r}$ = internal coordinate (transferred particle ↔ core of projectile)
- $\vec{R}$ = internal coordinate (transferred particle ↔ target)

These are related by:

$$\vec{r}_x = S_1\vec{r}_\alpha + T_1\vec{r}_\beta \quad \text{(transferred particle relative to target)}$$
$$\vec{r}_p = S_2\vec{r}_\alpha + T_2\vec{r}_\beta \quad \text{(transferred particle relative to ejectile)}$$

For stripping (ISTRIP=1), coefficients from source.mor lines 15876–15882.

### 6.2 Magnitudes

$$r_x = \sqrt{S_{1}^{2} \, r_{\alpha}^{2} + T_{1}^{2} \, r_{\beta}^{2} + 2 S_1 T_1 r_{\alpha} r_{\beta} \cos\phi_{ab}}$$
$$r_p = \sqrt{S_{2}^{2} \, r_{\alpha}^{2} + T_{2}^{2} \, r_{\beta}^{2} + 2 S_2 T_2 r_{\alpha} r_{\beta} \cos\phi_{ab}}$$

### 6.3 (U,V) Coordinate Transformation

$$U = \frac{r_\alpha + r_\beta}{2}, \quad V = r_\alpha - r_\beta, \quad \text{Jacobian} = 1$$

Domain: $U \in [0, U_{max}]$, $V \in [-2U, +2U]$ (ensures $r_\alpha > 0$, $r_\beta > 0$).

---

## 7. Integration Grid (GRDSET)

### 7.1 CUBMAP Adaptive Quadrature

Ptolemy uses cubic/rational-sinh mapped Gauss-Legendre points to concentrate quadrature near the integrand peak.

**Rational-sinh formula (MAPTYP=2):** Given GL nodes $t_k \in [-1,1]$ with weights $w_k$:

$$x_k = \frac{-A + (C/\gamma)\sinh(\tau t_k)}{B - (D/\gamma)\sinh(\tau t_k)}, \quad \tau = \text{arcsinh}(\gamma)$$

Points pile up near $x_{mid}$.

**Default parameters:**
| Axis | Map type | $\gamma$ | Midpoint |
|---|---|---|---|
| Sum U=(ri+ro)/2 | MAPTYP=2 (rational-sinh) | 2.0 | SUMMID (adaptive) |
| Dif V=ri-ro | MAPTYP=1 (cubic-sinh) | 12.0 | 0 (symmetric) |
| Phi φ_ab | MAPTYP=2 (rational-sinh) | ~0 (≈linear) | PHI0/2 |

### 7.2 Adaptive Phi Cutoff

For each (ri,ro) pair, Ptolemy does a 2-pass scan:

**Pass 1:** Scan $\cos\phi$ from +1 downward. Find last index where $|\text{BSPROD}| > \text{RVRLIM}/(r_i \cdot r_o)$.

**Pass 2:** Set $\phi_0 = \arccos(1 - 2(\text{IEND}-1)^2/N_{scan}^2)$, then integrate $\phi \in [0, \phi_0]$ with NPPHI GL points.

---

## 8. Radial Integrals (INELDC)

### 8.1 Transfer Integral

$$I_{L_i, L_o, L_x}^{J_{\pi}, J_{\pi}'} = \int_0^\infty dr_i \int_0^\infty dr_o \int_{-1}^{1} d(\cos\phi_{ab}) \; \chi_{L_i}^{J_{\pi}}(r_i) \cdot \mathcal{K}(r_i, r_o, \phi_{ab}) \cdot \chi_{L_o}^{J_{\pi}' \ast}(r_o)$$

where the kernel is:

$$\mathcal{K} = \phi_T(r_x) \cdot V(r_p) \cdot \phi_P(r_p) \cdot \mathcal{A}_{12}(r_i, r_o, \phi_{ab})$$

### 8.2 Angular Coupling Kernel (A12)

$$\mathcal{A}_{12}(r_i, r_o, \phi_{ab}) = \sum_{M_T, M_U} C_{M_T, M_U} \cos(M_T \phi_T - M_U \phi_{ab})$$

where $C_{M_T, M_U}$ involves Wigner d-matrix elements (xlam), 3j-symbols, and $\sqrt{2L_o+1}$ normalization factors.

### 8.3 INELDC Loop Structure

```
for LI (incoming L):
  compute χ_Li for incoming channel
  for LO (outgoing L):
    compute χ_Lo for outgoing channel
    for Lx (angular momentum transfer):
      compute CG angular coupling (A12)
      accumulate radial integral I(Li, Lo, Lx)
call SFROMI to convert I → S-matrix element
```

---

## 9. Transfer S-Matrix (SFROMI)

### 9.1 Assembly

$$S_{\text{SFROMI}}(L_i, L_o, L_x, J_{\pi}, J_{\pi}') = \text{FACTOR} \cdot \text{ATERM} \cdot \frac{i^{L_i+L_o+2L_x+1}}{\sqrt{2L_i+1}} \cdot I_{L_i, L_o, L_x}^{J_{\pi}, J_{\pi}'}$$

### 9.2 FACTOR

$$\text{FACTOR} = 2\sqrt{\frac{k_a k_b}{E_{\text{cm}}^{a} \, E_{\text{cm}}^{b}}}$$

### 9.3 ATERM (Spectroscopic Amplitude)

$$\text{ATERM}(L_x) = \sqrt{\frac{2j_B+1}{2j_A+1}} \cdot \sqrt{2L_x+1} \cdot \mathcal{S}_{\text{proj}} \cdot \mathcal{S}_{\text{target}} \cdot W(j_B, l_T, j_A, l_P; j_n, L_x)$$

### 9.4 Phase Sign

If $(L_i + L_o + 2L_x + 1) \bmod 4 \geq 2$: flip sign.
If $(L_i + L_o + 2L_x + 1) \bmod 2 = 1$: multiply by $i$.

### 9.5 9-J Structure

The full SFROMI includes a double 9-J symbol loop coupling spin quantum numbers of incoming/outgoing channel J-values:

$$S_{\text{SFROMI}} = \text{FACTOR} \cdot \mathcal{N}_{9J} \cdot \text{SAV9J} \cdot \text{TEMP} \cdot I$$

---

## 10. Cross Section

### 10.1 BETCAL — Partial Wave Amplitudes

$$\beta(L_o, L_x, M_x) = \frac{1}{2k_a} \sum_{L_i} (2L_i+1) \cdot C(L_i, 0, L_x, M_x; L_o, M_x) \cdot e^{i(\sigma_{L_i}^{(a)} + \sigma_{L_o}^{(b)})} \cdot S_{\text{SFROMI}}(L_i, L_o, L_x)$$

### 10.2 AMPCAL — Scattering Amplitude

$$F^{(M_x)}(\theta) = \sum_{L_o} \beta(L_o, L_x, M_x) \cdot P_{L_o}^{M_x}(\cos\theta)$$

### 10.3 XSECTN — Differential Cross Section

$$\frac{d\sigma}{d\Omega}(\theta) = 10 \cdot \sum_{L_x, M_x} f_{M_x} \cdot \left|F^{(M_x)}(\theta)\right|^{2}$$

where $f_{M_x} = 1$ for $M_x=0$, $f_{M_x} = 2$ for $M_x > 0$ (sum over $\pm M_x$), and the factor 10 converts fm² to mb.

---

## 11. Comparison: Handbook vs Ptolemy

| Aspect | Handbook | Ptolemy |
|---|---|---|
| **ODE** | Eq. 5: standard Schrödinger | Same in WAVELJ |
| **Regular BC** | $u \sim r^{l+1}$ | `u[1] = h^{L+1}`, Numerov forward |
| **Asymptotic BC** | Eq. 27: $\chi \to (e^{i\sigma}/2i)(S H^{+} - H^{-})$ | Same; ALPHAR/ALPHAI enforce this |
| **S-matrix extraction** | Wronskian | Two-point Wronskian at $r_N$, $r_{N-1}$ |
| **Coulomb phase** | $\arg\Gamma(L+1+iη)$ | SIGZRO array |
| **Coulomb amplitude** | §2.3, Eq. 10 | AMPCAL, exact formula |
| **Nuclear amplitude** | $f_N = \frac{1}{2ik}\sum(2L+1)e^{2i\sigma}(S-1)P_L$ | BETCAL + AMPCAL |
| **Approximation** | Both ZR and FR | Finite-range only |
| **T-matrix form** | Eq. 33 | Same. Prior/Post form |
| **Coordinate transforms** | Fig. 9 with S1, T1, S2, T2 | GRDSET lines 15876–15882 |
| **Radial integral** | Eq. 5.8: $\iint \chi_{\beta}^{\ast} F \chi_{\alpha} \, dr_{\alpha} \, dr_{\beta}$ | INELDC: 2D Gauss quadrature |
| **Kinematic factor** | $\mu_{\alpha}\mu_{\beta}/\hbar^4 \cdot k_{\beta}/k_{\alpha}$ | SFROMI: $2\sqrt{k_{\alpha} k_{\beta} / (E_{\text{cm}}^{(\alpha)} E_{\text{cm}}^{(\beta)})}$ (equivalent) |
| **Angular coupling** | 9j-symbol + CG | A12 subroutine |
| **Phase convention** | $i^{L_{\alpha} - L_{\beta} - l}$ | $i^{L_a + L_o + 2 L_x + 1}$ (equivalent) |
| **Spin-orbit** | §2.9: $\Lambda(L,s,J)$ | Same formula |
| **Units** | fm² | ×10 for mb/sr |

---

*Generated by Dudu 🐱 — sources: source.mor, handbook*
