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
8. [Transfer Cross Section: Overview](#8-transfer-cross-section-overview) ← **start here for the big picture**
9. [Radial Integrals (INELDC)](#9-radial-integrals-ineldc)
10. [Transfer S-Matrix (SFROMI)](#10-transfer-s-matrix-sfromi)
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

$$S_L = \frac{u_{N-1} \cdot H^{+}_{L}(r_N) - u_N \cdot H^{+}_{L}(r_{N-1})}{u_{N-1} \cdot H^{-}_{L}(r_N) - u_N \cdot H^{-}_L(r_{N-1})}$$

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
$$\text{Diff} = \left.\frac{u'}{u}\right|_\text{out} - \left.\frac{u'}{u}\right|_\text{in} \to 0$$

### 3.3 Depth Search Iteration

Iterate on potential depth $V$ to drive $\text{Diff} \to 0$:

1. **Node check**: If nodes > n, V too deep (decrease); if nodes < n, V too shallow (increase)
2. **Newton-Raphson**: Once node count is correct, update $V_{new} = V_{old} + \text{factor} \times \text{Diff}$
3. **Adaptive factor**: Starts at 2.0; halved when oscillation detected (sign of Diff flips)

---

## 4. DWBA Transfer Reaction

### 4.1 The Transition Amplitude

For a transfer reaction $a (= b + x) + A \to b + B (= A + x)$:

$$T_{ba} = \int d^3 r_b \int d^3 r_a \, \chi_{b}^{(-)\ast}(\mathbf{r}_{b}) \langle \psi_B \psi_b | V | \psi_A \psi_a \rangle \chi_{a}^{(+)}(\mathbf{r}_a)$$

Separating internal states from radial motion:

$$T_{ba} \propto \int d^3 r_b \int d^3 r_a \, \chi_{b}^{\ast}(\mathbf{r}_{b}) \, \phi_{Bx}^{\ast}(\mathbf{r}_{x}) \, V(\mathbf{r}_{bx}) \, \phi_{ax}(\mathbf{r}_{bx}) \, \chi_{a}(\mathbf{r}_a)$$

- $\phi_{Bx}$: Bound state wavefunction of x in B (target bound state)
- $\phi_{ax}$: Bound state wavefunction of x in a (projectile bound state)

### 4.2 Partial Wave Expansion

Expanding all wavefunctions in partial waves:

$$T \sim \sum_{L_a, L_b} \int r_{a}^{2} \, dr_a \int r_b^{2} \, dr_b \, \chi_{L_b}(r_b) \, K_{L_a L_b}(r_a, r_b) \, \chi_{L_a}(r_a)$$

The kernel $K$ involves angular integration:

$$K(r_a, r_b) = \int_{-1}^{1} d(\cos \theta) \, \phi_{Bx}^{\ast}(r_x) \, V(r_{bx}) \, \phi_{ax}(r_{bx}) \, P_L(\cos \theta)$$

where $r_x$ and $r_{bx}$ depend on $r_a, r_b, \theta$ via the law of cosines.

### 4.3 Stripping vs Pickup

Ptolemy handles both directions of transfer using the **same integral code** with different coordinate mappings.

| ISTRIP | Direction | Example | Definition |
|---|---|---|---|
| +1 | Stripping | (d,p), (t,d), (³He,d) | Projectile loses particle: $a = b + x$ |
| −1 | Pickup | (p,d), (d,t), (d,³He) | Projectile gains particle: $b = a + x$ |
| 0 | Inelastic | (p,p'), (d,d') | No mass transfer |

**Detection is automatic:** Ptolemy compares $m_a$ vs $m_b$. If $m_a > m_b$ → stripping; if $m_a < m_b$ → pickup.

For **pickup** (ISTRIP=−1), Ptolemy flips the S1/T1/S2/T2 mapping (source.mor lines 15889–15896):

```
T2 = S2;  S2 = -S1;  S1 = T1;  T1 = -S2;  PHISGN = -1
```

It also adjusts:
- Q-value sign: `Q = -Q` (line 4490)
- ATERM stat factor: `(JB+1)/(JA+1)` instead of `(JBIGB+1)/(JBIGA+1)` (line 25631)
- Bound state assignments: which wavefunction is "projectile" vs "target" (lines 27893–27901)

The integral machinery (INELDC, SFROMI, BETCAL, XSECTN) is identical — all the stripping/pickup physics is encoded in S1/T1/S2/T2 and the bound state assignments.

### 4.4 Multi-Nucleon Transfer

Ptolemy treats the transferred cluster $x$ as a **single entity**, regardless of how many nucleons it contains. For example:

| Reaction | x (cluster) | $m_x$ | $j_x$ | Notes |
|---|---|---|---|---|
| (d,p) | n | 1 | 1/2 | Single neutron |
| (t,p) | 2n | 2 | 0 | Dineutron (J=0 assumed) |
| (³He,p) | 2n (or d) | 2 | 0 or 1 | Depends on model |
| (³He,d) | p | 1 | 1/2 | Single proton |
| (α,d) | 2n | 2 | 0 | Dineutron cluster |
| (t,d) | n | 1 | 1/2 | Single neutron pickup |

**For dineutron/diproton (A=2, Z=0 or Z=2):** Ptolemy automatically assigns $J_x = 0$ (source.mor line 7064–7068). This is the standard assumption — the two transferred nucleons are in a relative $^1S_0$ state.

**What this means physically:** The (t,p) reaction is treated as transferring a di-neutron cluster with:
- **Projectile bound state:** 2n inside the triton ($t = p + 2n$)
- **Target bound state:** 2n inside the residual nucleus ($B = A + 2n$)

Both bound states get their own Woods-Saxon potential, orbital quantum numbers, and binding energies. The DWBA integral is then the same as for single-nucleon transfer — Ptolemy does not do sequential (two-step) transfer.

> **Limitation:** This "cluster DWBA" approach assumes simultaneous transfer. For reactions where sequential transfer dominates (e.g., some (t,p) reactions at higher energies), coupled-channels or second-order DWBA would be needed. Ptolemy's coupled-channels capability is limited and rarely used for this purpose.

---

## 5. Zero-Range vs Finite-Range

### Zero-Range Approximation

Assumes the interaction and projectile bound state are pointlike:

$$V(\mathbf{r}_{bx}) \, \phi_{ax}(\mathbf{r}_{bx}) \approx D_0 \, \delta(\mathbf{r}_{bx})$$

This collapses the 6D integral to a 3D integral:

$$T_{\text{ZR}} \propto \int d^3 r \, \chi_{b}^{\ast}(\beta \mathbf{r}) \, \phi_{Bx}^{\ast}(\mathbf{r}) \, \chi_a(\mathbf{r})$$

Fast but inaccurate for heavy ion reactions, high energies, and states where finite size matters. See §9.4 for how the finite-range integral collapses to ZR.

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

### 6.1 The Three-Body Problem

In a transfer reaction $A(a,b)B$ where $a = b + x$ and $B = A + x$, there are four particles but only two independent relative coordinates. The challenge is that the **bound state** coordinates (x relative to its core) differ from the **scattering** coordinates (projectile relative to target).

We define:
- $\mathbf{r}_\alpha$ = relative coordinate between a and A (entrance channel)
- $\mathbf{r}_\beta$ = relative coordinate between b and B (exit channel)
- $\mathbf{r}_T$ = position of x relative to A (target bound state coordinate)
- $\mathbf{r}_P$ = position of x relative to b (projectile bound state coordinate)

The key insight: the bound state coordinates are **linear combinations** of the scattering coordinates:

$$\mathbf{r}_T = S_1 \mathbf{r}_\alpha + T_1 \mathbf{r}_\beta$$
$$\mathbf{r}_P = S_2 \mathbf{r}_\alpha + T_2 \mathbf{r}_\beta$$

### 6.2 Derivation of S1, T1, S2, T2

For **stripping** ($a \to b + x$, ISTRIP=1), the coefficients follow from center-of-mass relations. Define mass ratios:

$$\rho_1 = \frac{m_x}{m_b}, \qquad \rho_2 = \frac{m_x}{m_A}$$

Then (Fortran source.mor lines 15876–15881):

$$\text{TEMP} = \frac{1}{\rho_1 + \rho_2(1 + \rho_1)}$$

$$S_1 = (1 + \rho_1)(1 + \rho_2) \cdot \text{TEMP}$$

$$T_1 = -(1 + \rho_2) \cdot \text{TEMP}$$

$$S_2 = (1 + \rho_1) \cdot \text{TEMP}$$

$$T_2 = -S_1$$

**Physical meaning:**
- $S_1$ and $T_1$ project channel coordinates onto the **target bound state** coordinate $r_T$ (x inside B)
- $S_2$ and $T_2$ project onto the **projectile bound state** coordinate $r_P$ (x inside a)
- $T_2 = -S_1$ is an exact identity from momentum conservation

**Example:** For 16O(d,p)17O with $m_x = m_n$, $m_b = m_p$, $m_A = m_{16\text{O}}$:

$$\rho_1 = \frac{m_n}{m_p} \approx 1.001, \quad \rho_2 = \frac{m_n}{m_{16\text{O}}} \approx 0.063$$

$$S_1 \approx 1.887, \quad T_1 \approx -0.944, \quad S_2 \approx 0.943, \quad T_2 \approx -1.887$$

For **pickup** (ISTRIP=−1), the roles of a and B are swapped, and the coefficients are rearranged (source.mor lines 15889–15896).

### 6.3 Magnitudes

Since the bound state coordinates are vector sums of scattering coordinates, their magnitudes involve the angle $\phi_{ab}$ between the two channel vectors:

$$r_T = \sqrt{S_{1}^{2} \, r_{\alpha}^{2} + T_{1}^{2} \, r_{\beta}^{2} + 2 S_1 T_1 \, r_\alpha \, r_\beta \cos\phi_{ab}}$$

$$r_P = \sqrt{S_{2}^{2} \, r_{\alpha}^{2} + T_{2}^{2} \, r_{\beta}^{2} + 2 S_2 T_2 \, r_\alpha \, r_\beta \cos\phi_{ab}}$$

Note that $S_1 T_1 < 0$ and $S_2 T_2 < 0$, so the cross-term reduces the magnitude — the bound state coordinates are always shorter than the scattering coordinates.

### 6.4 (U,V) Coordinate Transformation

For numerical integration, Ptolemy uses sum/difference coordinates:

$$U = \frac{r_{\alpha} + r_{\beta}}{2}, \quad V = r_\alpha - r_\beta, \quad \text{Jacobian} = 1$$

Domain: $U \in [0, U_{\max}]$, $V \in [-2U, +2U]$ (ensures $r_\alpha > 0$, $r_\beta > 0$).

This is advantageous because the integrand peaks near $r_\alpha \approx r_\beta$ (i.e., $V \approx 0$), allowing efficient quadrature concentration.

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

## 8. Transfer Cross Section: Overview

This section shows how the full DWBA transition amplitude factorizes into the components computed by Ptolemy's subroutines. The subsequent chapters (§9–§12) detail each piece.

### 8.1 The Full T-Matrix

Starting from the DWBA transition amplitude (§4.1):

$$T_{ba} = \int d^3 r_\beta \int d^3 r_\alpha \; \chi_b^{(-)\ast}(\mathbf{r}_\beta) \; \phi_T^\ast(\mathbf{r}_T) \; V(\mathbf{r}_P) \; \phi_P(\mathbf{r}_P) \; \chi_a^{(+)}(\mathbf{r}_\alpha)$$

This is a 6D integral involving:
- **Distorted waves** $\chi_a$, $\chi_b$ — scattering states in entrance/exit channels (§2)
- **Bound states** $\phi_T$, $\phi_P$ — bound wavefunctions for target and projectile vertices (§3)
- **Interaction** $V(\mathbf{r}_P)$ — the binding potential at the projectile vertex
- **Coordinate mapping** from scattering coordinates to bound state coordinates (§6)

### 8.2 From T-Matrix to Cross Section

Expand all wavefunctions in partial waves. The angular integrals produce Clebsch-Gordan coefficients, Racah W-coefficients, and 9-j symbols. The cross section factorizes into three stages:

**Step 1 — SFROMI:** Combine the radial integral with kinematic and spectroscopic factors to get a transfer S-matrix element:

$$S(L_i, L_o, L_x) = \text{FACTOR} \cdot \text{ATERM}(L_x) \cdot \frac{i^{L_i+L_o+2L_x+1}}{\sqrt{2L_i+1}} \cdot I_{L_i,L_o,L_x}^{J_\pi,J_\pi'}$$

where:
- FACTOR $= 2\sqrt{k_a k_b / (E_\text{cm}^a E_\text{cm}^b)}$ — kinematic factor (§10.2)
- ATERM — spectroscopic amplitude with Racah coefficient (§10.3)
- $I$ — the 3D radial-angular integral computed by INELDC (§9)

**Step 2 — BETCAL:** Sum over incoming partial waves $L_i$ with Clebsch-Gordan coupling and Coulomb phases:

$$\beta(L_o, L_x, M_x) = \frac{1}{2k_a} \sum_{L_i} (2L_i+1) \cdot C(L_i, 0; L_x, M_x | L_o, M_x) \cdot e^{i(\sigma_{L_i} + \sigma_{L_o})} \cdot S(L_i, L_o, L_x)$$

where $C(L_i, 0; L_x, M_x | L_o, M_x)$ is a Clebsch-Gordan coefficient coupling the incoming partial wave $L_i$ (with $M=0$) and the transferred angular momentum $(L_x, M_x)$ to the outgoing partial wave $L_o$.

**Step 3 — AMPCAL:** Sum over outgoing partial waves to get the scattering amplitude at angle $\theta$:

$$F^{(M_x)}(\theta) = \sum_{L_o} \beta(L_o, L_x, M_x) \cdot P_{L_o}^{M_x}(\cos\theta)$$

where $P_{L_o}^{M_x}$ is an associated Legendre polynomial.

**Step 4 — XSECTN:** Sum incoherently over the transferred angular momentum projection $M_x$:

$$\frac{d\sigma}{d\Omega}(\theta) = 10 \cdot \sum_{L_x, M_x} f_{M_x} \cdot \left|F^{(M_x)}(\theta)\right|^{2}$$

where:
- The factor **10** converts fm² → mb/sr
- $f_{M_x} = 1$ for $M_x = 0$, and $f_{M_x} = 2$ for $M_x > 0$ (because $M_x$ and $-M_x$ give identical contributions, so we sum only $M_x \geq 0$ and double)

### 8.3 The Factorization Map

The cross section computation splits into five independent stages, each handled by a dedicated subroutine:

```
Stage 1: BOUND (§3)              → φ_T(r), φ_P(r)           [bound state wavefunctions]
Stage 2: WAVELJ (§2)             → χ_Li(r), χ_Lo(r)         [distorted waves]
Stage 3: GRDSET + INELDC (§7,§9) → I(Li, Lo, Lx, Jπ, Jπ')  [radial integral]
Stage 4: SFROMI (§10)            → S(Li, Lo, Lx)            [transfer S-matrix]
Stage 5: BETCAL + AMPCAL (§8.2)  → F(θ), dσ/dΩ             [cross section]
```

**Why this factorization matters:**
- **Stages 1–2** are 1D problems (single ODE each) — cheap
- **Stage 3** is the expensive part — a 3D integral over $(r_\alpha, r_\beta, \phi_{ab})$ for each $(L_i, L_o, L_x, J_\pi)$ combination
- **Stages 4–5** are algebraic — just sums over angular momentum quantum numbers

### 8.4 Quantum Number Loops

The transfer cross section involves nested sums over angular momentum quantum numbers. Here is the complete loop structure:

**Outer (XSECTN):** Sum over $L_x$ (transferred angular momentum) and $M_x$ (projection)

**Middle (BETCAL):** For each $(L_x, M_x)$, sum over $L_i$ (incoming partial wave), with $L_o$ fixed by CG selection rule

**Inner (SFROMI):** For each $(L_i, L_o, L_x)$, sum over $(J_\pi, J_\pi')$ via 9-j symbols coupling channel spins

**Innermost (INELDC):** For each $(L_i, L_o, L_x, J_\pi, J_\pi')$, compute the 3D radial-angular integral

The total number of integrals scales as $N_L^2 \times N_{L_x} \times N_J^2$, where $N_L \sim 30\text{--}40$, $N_{L_x} \sim 1\text{--}5$, and $N_J \sim 2$ (for spin-1/2 channels).

---

## 9. Radial Integrals (INELDC)

### 9.1 From 6D to the Radial-Angular Form

The DWBA transition amplitude (§4.1) is a 6-dimensional integral over the entrance and exit channel coordinates. By expanding all wavefunctions in partial waves (spherical harmonics), the angular integrals can be done analytically, reducing the problem to radial integrals.

**Step 1:** Expand distorted waves in partial waves:

$$\chi_a(\mathbf{r}_{\alpha}) = \sum_{L_i, M_i} \frac{u_{L_i}(r_{\alpha})}{r_{\alpha}} Y_{L_i}^{M_i}(\hat{r}_\alpha), \qquad \chi_b(\mathbf{r}_{\beta}) = \sum_{L_o, M_o} \frac{u_{L_o}(r_{\beta})}{r_{\beta}} Y_{L_o}^{M_o}(\hat{r}_\beta)$$

**Step 2:** Expand bound states in spherical harmonics:

$$\phi_T(r_T) Y_{l_T}^{m_T}(\hat{r}_T), \qquad \phi_P(r_P) Y_{l_P}^{m_P}(\hat{r}_P)$$

**Step 3:** The angular integrals over $\hat{r}_\alpha$ and $\hat{r}_\beta$ produce Clebsch-Gordan coefficients and 3j-symbols, coupling $(L_i, l_T, l_P, L_o)$ to a transferred angular momentum $L_x$. Selection rules enforce $|L_i - L_o| \leq L_x \leq L_i + L_o$ and $|l_T - l_P| \leq L_x \leq l_T + l_P$.

**Step 4:** After angular reduction, the remaining integral is over three scalar variables $(r_\alpha, r_\beta, \phi_{ab})$:

$$I_{L_i, L_o, L_x}^{J_{\pi}, J_{\pi}'} = \int_0^{\infty} dr_{\alpha} \int_0^{\infty} dr_{\beta} \int_{-1}^{1} d(\cos\phi_{ab}) \; \chi_{L_i}^{J_{\pi}}(r_{\alpha}) \cdot \mathcal{K}(r_\alpha, r_{\beta}, \phi_{ab}) \cdot \chi_{L_o}^{J_{\pi}' \ast}(r_\beta)$$

The angle $\phi_{ab}$ is the angle between the entrance and exit channel coordinate vectors — it survives because the bound state coordinates depend on it through the law of cosines (§6.3).

### 9.2 The Kernel

$$\mathcal{K}(r_{\alpha}, r_{\beta}, \phi_{ab}) = \phi_T(r_T) \cdot V(r_P) \cdot \phi_P(r_P) \cdot \mathcal{A}_{12}(r_\alpha, r_\beta, \phi_{ab})$$

where:
- $\phi_T(r_T)$ = target bound state wavefunction, evaluated at $r_T(r_\alpha, r_\beta, \phi_{ab})$
- $V(r_P) \cdot \phi_P(r_P)$ = interaction × projectile bound state, evaluated at $r_P(r_\alpha, r_\beta, \phi_{ab})$
- $\mathcal{A}_{12}$ = angular coupling coefficient from the partial wave reduction

The kernel is **sharply peaked** near $\cos\phi_{ab} \to +1$ (collinear geometry), because this is where $r_P \to 0$ (the transferred particle is closest to the projectile core), making $V(r_P) \cdot \phi_P(r_P)$ large.

### 9.3 Angular Coupling Kernel (A12)

$$\mathcal{A}_{12}(r_\alpha, r_\beta, \phi_{ab}) = \sum_{M_T, M_U} C_{M_T, M_U} \cos(M_T \phi_T - M_U \phi_{ab})$$

where $C_{M_T, M_U}$ involves Wigner d-matrix elements (xlam), 3j-symbols, and $\sqrt{2L_o+1}$ normalization factors. This encodes the geometric coupling between the partial waves of the scattering and bound states.

### 9.4 Connection to Zero-Range

In the **zero-range limit**, the projectile bound state and interaction collapse to a delta function:

$$V(r_P) \cdot \phi_P(r_P) \to D_0 \, \delta^{(3)}(\mathbf{r}_P)$$

This forces the projectile coordinate to zero. From the mapping $\mathbf{r}_P = S_2 \mathbf{r}_\alpha + T_2 \mathbf{r}_\beta$, this implies:

$$\mathbf{r}_\beta = -\frac{S_2}{T_2} \mathbf{r}_{\alpha} = \frac{S_2}{S_1} \mathbf{r}_\alpha \quad (\text{since } T_2 = -S_1)$$

So the exit channel coordinate is locked proportional to the entrance channel coordinate. The $\phi_{ab}$ integral collapses (the delta function fixes the angle), and the target coordinate becomes a function of a single radial variable:

$$r_T = |S_1 - T_1 S_2/S_1| \cdot r_{\alpha} = \left|S_1 + T_1 \frac{S_2}{S_1}\right| r_\alpha$$

The 6D integral reduces to a **single radial integral**:

$$T_{\text{ZR}} \propto D_0 \int_0^{\infty} dr \; \chi_b^{\ast}\!\left(\frac{S_2}{S_1} r\right) \phi_T(r) \, \chi_a(r)$$

This is computationally trivial but ignores the finite spatial extent of the projectile — leading to 10–30% errors for (d,p) reactions where the deuteron has a large radius (~4 fm).

### 9.5 INELDC Loop Structure

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

## 10. Transfer S-Matrix (SFROMI)

### 10.1 Assembly

$$S_{\text{SFROMI}}(L_i, L_o, L_x, J_{\pi}, J_{\pi}') = \text{FACTOR} \cdot \text{ATERM} \cdot \frac{i^{L_i+L_o+2L_x+1}}{\sqrt{2L_i+1}} \cdot I_{L_i, L_o, L_x}^{J_{\pi}, J_{\pi}'}$$

### 10.2 FACTOR

$$\text{FACTOR} = 2\sqrt{\frac{k_a k_b}{E_{\text{cm}}^{a} \, E_{\text{cm}}^{b}}}$$

### 10.3 ATERM (Spectroscopic Amplitude)

**Fortran BSSET (source.mor line 25634):**

$$\text{ATERM}(L_x) = \sqrt{\frac{J_B'+1}{J_A'+1}} \cdot \sqrt{2L_x+1} \cdot \mathcal{S}_{\text{proj}} \cdot \mathcal{S}_{\text{target}} \cdot W(l_T, j_T, l_P, j_P; j_x, L_x)$$

where $J_A' = 2J_A$ and $J_B' = 2J_B$ are the doubled nuclear spin quantum numbers (Ptolemy convention), and $W$ is the Racah coefficient (related to 6-j by a phase).

**Sign:** $(-1)^{(j_x - j_P + 2(l_P + l_T))/2 + 1}$ flips ATERM when this exponent is odd.

> **✅ C++ code note (no double-count):**
> The faithful path (`grdset_ineldc_faithful.cpp`) applies FACTOR, ATERM (with Racah), phase, and `1/√(2Li+1)` when building `TransferSMatrix` entries. In `xsectn.cpp`, there are two loops over `TransferSMatrix`: the first (line ~133) redundantly re-applies FACTOR/ATERM, but its results are **overwritten** by the second loop (line ~255), which correctly uses `elem.S` directly (already containing FACTOR×ATERM) and feeds it into the 9-J accumulation without re-applying these factors.
>
> **Verified:** Comparing C++ vs Fortran DCS for both 16O(d,p)17O and 206Hg(d,p)207Hg shows an angle-dependent error of ~8-10%, confirming the error is in the radial integral (GRDSET/INELDC), not in ATERM/FACTOR. A double-counted ATERM would produce a constant scale factor independent of angle.

### 10.4 Phase Sign

If $(L_i + L_o + 2L_x + 1) \bmod 4 \geq 2$: flip sign.
If $(L_i + L_o + 2L_x + 1) \bmod 2 = 1$: multiply by $i$.

### 10.5 9-J Structure

When spin-orbit potentials are present, the transfer S-matrix involves a double sum over channel spin states via 9-J symbols. For each radial integral element with quantum numbers $(L_i, L_o, L_x^{(\text{phys})}, J_\pi, J_\pi')$, the accumulation into the cross-section S-matrix slot $(L_i, L_o, L_x, J_P)$ proceeds as:

$$S(k, L_i) \mathrel{+}= \sum_{j_{\pi i}, j_{\pi o}} \text{9J}_1 \cdot (SR + i\,SI) \cdot \text{9J}_2 \cdot (-1)^{L_x + L_x^{(\text{phys})}}$$

where:

**First 9-J symbol** (SAV9J in the code) couples the physical quantum numbers to the channel spin:

$$\text{9J}_1 = \sqrt{(2j_{\pi i}+1)(2j_{\pi o}+1)(2L_x^{(\text{phys})}+1)(J_{BP}+1)} \;\left\{ \begin{matrix} J_{BT}/2 & L_x^{(\text{phys})} & J_{BP}/2 \\ j_{\pi i}/2 & L_i & J_A/2 \\ j_{\pi o}/2 & L_o & J_B/2 \end{matrix} \right\}$$

**Second 9-J symbol** (TEMP2) recouples into the cross-section angular momentum $(L_x, J_P)$:

$$\text{9J}_2 = \sqrt{(2j_{\pi i}+1)(2j_{\pi o}+1)(2L_x+1)(J_P+1)} \;\left\{ \begin{matrix} J_{BT}/2 & L_x & J_P/2 \\ j_{\pi i}/2 & L_i & J_A/2 \\ j_{\pi o}/2 & L_o & J_B/2 \end{matrix} \right\}$$

Here:
- $j_{\pi i}$, $j_{\pi o}$ are the channel spin couplings ($J_\pi = L \otimes S_{\text{channel}}$) in the incoming/outgoing channels
- $J_A$, $J_B$ are the doubled target/residual spins
- $J_{BT}$, $J_{BP}$ are the doubled transferred particle total angular momenta
- $(SR, SI)$ is the complex S-matrix element from GRDSET/INELDC (with FACTOR, ATERM, and phase already applied)
- The phase $(-1)^{L_x + L_x^{(\text{phys})}}$ accounts for the parity difference when $L_x \neq L_x^{(\text{phys})}$

When the two 9-J symbols share the same arguments (i.e., $L_x = L_x^{(\text{phys})}$ and $J_P = J_{BP}$), the second 9-J equals the first (optimization in the code).

---

## 11. Comparison: Handbook vs Ptolemy

| Aspect | Handbook | Ptolemy |
|---|---|---|
| **ODE** | §2.1, Eq. 2.5 | Same in WAVELJ |
| **Regular BC** | $u \sim r^{l+1}$ | `u[1] = h^{L+1}`, Numerov forward |
| **Asymptotic BC** | §2.4, Eq. 2.27 | Same; ALPHAR/ALPHAI enforce this |
| **S-matrix extraction** | §2.7: Wronskian | Two-point Wronskian at $r_N$, $r_{N-1}$ |
| **Coulomb phase** | $\arg\Gamma(L+1+i\eta)$ | SIGZRO array |
| **Coulomb amplitude** | §2.3, Eq. 2.10 | AMPCAL, exact formula |
| **Nuclear amplitude** | §2.4: partial wave sum | BETCAL + AMPCAL |
| **Approximation** | Both ZR and FR | Finite-range only |
| **T-matrix form** | §5.2, Eq. 5.33 | Same. Prior/Post form |
| **Coordinate transforms** | §5.1, Fig. 5.1 | GRDSET lines 15876–15882 |
| **Radial integral** | §5.8, Eq. 5.34 | INELDC: 2D Gauss quadrature |
| **Kinematic factor** | §5.5, Eq. 5.36 | SFROMI: equivalent form |
| **Angular coupling** | §5.3: 9j-symbol + CG | A12 subroutine |
| **Phase convention** | $i^{L_\alpha - L_\beta - l}$ | $i^{L_a + L_o + 2 L_x + 1}$ (equivalent) |
| **Spin-orbit** | §2.9, Eq. 2.18 | Same formula |
| **Units** | fm² | ×10 for mb/sr |

---

*Generated by Dudu 🐱 — sources: source.mor, handbook*
