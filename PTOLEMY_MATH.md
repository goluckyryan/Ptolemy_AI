# Ptolemy DWBA: Complete Mathematical Reference

**Reaction:** 33Si(d,p)34Si at E_lab = 20 MeV, L=2 transfer (0d3/2 in 34Si)
**Last updated:** 2026-03-14

---

## 1. Kinematics

The reaction is A(a,b)B where:
- a = deuteron (d), b = proton (p)
- A = 33Si, B = 34Si
- Transferred nucleon: neutron n (lT=2, jT=3/2 in 34Si; lP=0, jP=1/2 in d)

### Center-of-mass energies
$$E_{cm}^a = E_{lab} \cdot \frac{M_A}{M_a + M_A} = 20 \cdot \frac{32.97}{1.01+32.97} = 18.843 \text{ MeV}$$

$$E_{cm}^b = E_{cm}^a + Q = 18.843 + 5.3241 = 24.167 \text{ MeV}$$

### Wave numbers
$$k_a = \sqrt{\frac{2\mu_a E_{cm}^a}{\hbar^2}}, \quad k_b = \sqrt{\frac{2\mu_b E_{cm}^b}{\hbar^2}}$$

Numerical: $k_a = 1.3109\ \text{fm}^{-1}$, $k_b = 1.0634\ \text{fm}^{-1}$

### Sommerfeld parameters
$$\eta_a = \frac{Z_a Z_A e^2 \mu_a}{\hbar^2 k_a} = 0.7043, \quad \eta_b = \frac{Z_b Z_B e^2 \mu_b}{\hbar^2 k_b} = 0.4516$$

---

## 2. Optical Model Potentials

### General Wood-Saxon form
$$U(r) = -V f(r, r_0, a) - i W_V f(r, r_i, a_i) + i W_S \frac{d}{dr}f(r, r_{si}, a_{si}) + V_{SO} \frac{1}{r}\frac{d}{dr}f(r, r_{so}, a_{so}) \, \vec{L}\cdot\vec{S} + V_C(r)$$

where $f(r, r_0, a) = \left[1 + \exp\!\left(\frac{r - r_0 A^{1/3}}{a}\right)\right]^{-1}$

### Coulomb potential
$$V_C(r) = \begin{cases} \frac{Z_p Z_T e^2}{2R_C}\left(3 - \frac{r^2}{R_C^2}\right) & r \leq R_C \\ \frac{Z_p Z_T e^2}{r} & r > R_C \end{cases}, \quad R_C = r_{C0} A^{1/3}$$

### Incoming channel (d + 33Si)
| Parameter | Value |
|---|---|
| V | 89.72 MeV |
| r0 | 1.1496 fm |
| a | 0.7594 fm |
| W_V | 2.348 MeV |
| r_i | 1.3361 fm |
| a_i | 0.5342 fm |
| W_S | 10.218 MeV |
| r_si | 1.3814 fm |
| a_si | 0.7299 fm |
| V_SO | 3.557 MeV |
| r_so | 0.972 fm |
| a_so | 1.011 fm |
| r_C | 1.303 fm |

### Outgoing channel (p + 34Si)
| Parameter | Value |
|---|---|
| V | 50.639 MeV |
| r0 | 1.179 fm |
| a | 0.673 fm |
| W_V | 2.397 MeV |
| r_i | 1.179 fm |
| a_i | 0.673 fm |
| W_S | 8.142 MeV |
| r_si | 1.291 fm |
| a_si | 0.536 fm |
| V_SO | 5.274 MeV |
| r_so | 0.986 fm |
| a_so | 0.590 fm |
| r_C | 1.301 fm |

---

## 3. Distorted Waves (WAVELJ)

Solve the radial Schrödinger equation with Numerov method:

$$\left[-\frac{\hbar^2}{2\mu}\frac{d^2}{dr^2} + \frac{\hbar^2}{2\mu}\frac{L(L+1)}{r^2} + U(r) - E_{cm}\right] u_L(r) = 0$$

### Effective potential (Numerov kernel)
$$f(r) = \frac{2\mu}{\hbar^2}\left[U(r) - E_{cm}\right] + \frac{L(L+1)}{r^2}$$

Note: nuclear potential **adds** to $k^2$ (attractive → negative contribution to $f$).

### Modified Numerov recursion
$$u_{i+1} = \frac{2\left(1 - \frac{5h^2}{12}f_i\right)u_i - \left(1 + \frac{h^2}{12}f_{i-1}\right)u_{i-1}}{1 + \frac{h^2}{12}f_{i+1}}$$

with $h=0.1$ fm, $u_0=0$, $u_1=h^{L+1}$.

### S-matrix extraction
At $r_{match}$ (beyond nuclear range), match to Coulomb asymptotic:
$$u_L(r) \to C\left[F_L(\eta, kr) - S_L \cdot G_L(\eta, kr)\right] \cdot e^{i\sigma_L}$$

$$S_L = \frac{u_{N-1} \cdot H^+_{L}(r_N) - u_N \cdot H^+_L(r_{N-1})}{u_{N-1} \cdot H^-_L(r_N) - u_N \cdot H^-_L(r_{N-1})}$$

where $H^\pm_L = G_L \pm iF_L$ and $F_L, G_L$ from RCWFN (Steed's method).

### Coulomb phases
$$\sigma_L = \arg\left[\Gamma(L+1+i\eta)\right] = \sum_{n=1}^{L} \arctan\!\left(\frac{\eta}{n}\right)$$

---

## 4. Bound State Wavefunctions (BOUND)

Same Numerov integration with real WS potential adjusted to give the correct binding energy $E_B$:

$V_{sol}$ found by bisection so that $u_L(r)$ has exactly $n$ nodes and $u_L \to e^{-\kappa r}$ for $r \to \infty$, where $\kappa = \sqrt{2\mu |E_B|/\hbar^2}$.

### Normalization
$$\int_0^\infty u_L^2(r)\, dr = 1, \quad \phi_L(r) = u_L(r)/r$$

### Target bound state (n in 34Si, 0d3/2)
- $E_B = -7.549$ MeV, $\kappa = 0.5939$ fm$^{-1}$, $l=2$, $j=3/2$
- $V_{sol} = 47.166$ MeV, $r_0=1.25$ fm, $a=0.65$ fm

### Projectile bound state (n in d, 0s1/2)
- $E_B = -2.225$ MeV, $\kappa = 0.2315$ fm$^{-1}$, $l=0$, $j=1/2$
- $V_{sol} = 62.942$ MeV, $r_0=1.25$ fm, $a=0.65$ fm
- AV18 spectroscopic amplitude: $\mathcal{S}_{AV18}^{1/2} = 0.97069$

---

## 5. Coordinate Geometry (BSPROD)

For the 3-body system (target core A, transferred neutron n, projectile b=p):

Let $\vec{r}_i$ = incoming Jacobi coordinate (d–33Si), $\vec{r}_o$ = outgoing Jacobi coordinate (p–34Si).

### Linear transformation coefficients

$$\vec{r}_x = S_1 \vec{r}_i + T_1 \vec{r}_o \quad \text{(n position relative to A)}$$
$$\vec{r}_p = S_2 \vec{r}_i + T_2 \vec{r}_o \quad \text{(n position relative to b=p)}$$

where (for stripping d→p+n):
$$S_1 = \frac{m_b}{m_b + m_n} = \frac{M_p}{M_p + M_n} = \frac{1.008}{2.016} \cdot 2 = 1.9429$$

$$T_1 = -\frac{M_B}{M_B + M_n} \cdot \frac{M_A}{M_a} \approx -1.8857$$

$$S_2 = \frac{M_A}{M_B} \cdot \frac{m_n}{m_b + m_n} \approx 0.9714, \quad T_2 = -S_1 = -1.9429$$

### Magnitudes
$$r_x = \sqrt{S_1^2 r_i^2 + T_1^2 r_o^2 + 2 S_1 T_1 r_i r_o \cos\phi_{ab}}$$
$$r_p = \sqrt{S_2^2 r_i^2 + T_2^2 r_o^2 + 2 S_2 T_2 r_i r_o \cos\phi_{ab}}$$

### Angle of rx in (ri, phi_ab) frame
$$\cos\phi_T = \frac{T_1 r_o + S_1 r_i \cos\phi_{ab}}{r_x}$$

### phi×V clipping (Ptolemy BSPROD ITYPE=2)
Ptolemy clips the **combined** product $\phi_T(r_x) \cdot V_{nA}(r_x)$ at its spatial maximum, not each factor separately:
$$[\phi_T \cdot V]_{clipped}(r) = \begin{cases} [\phi_T \cdot V](r) & r \leq r_{clip} \\ [\phi_T \cdot V]_{max} & r < r_{clip} \end{cases}$$

---

## 6. Angular Coupling Kernel (A12)

The angular coupling coefficient for the transfer matrix element:

$$\mathcal{A}_{12}(r_i, r_o, \phi_{ab}) = \sum_{M_T, M_U} C_{M_T, M_U} \cos(M_T \phi_T - M_U \phi_{ab})$$

where the coefficients are:
$$C_{M_T, M_U} = X_N \cdot \lambda(l_T, |M_T|) \cdot \lambda(l_P, 0) \cdot \begin{pmatrix} l_T & l_P & L_x \\ M_T & 0 & -M_T \end{pmatrix} \cdot \lambda(L_i, |M_U|) \cdot \lambda(L_o, |M_T - M_U|) \cdot \sqrt{2L_o+1}$$

$$\times \begin{pmatrix} L_i & L_o & L_x \\ M_U & M_T-M_U & -M_T \end{pmatrix} \cdot \epsilon(M_T) \cdot \epsilon(M_U)$$

with $X_N = \frac{1}{2}\sqrt{(2L_i+1)(2l_T+1)(2l_P+1)}$

### xlam (Wigner d-matrix)
$$\lambda(L, |M|) = d^L_{|M|,0}(\pi/2)$$

Computed via Ptolemy's alternating-M recurrence (NOT a simple closed form):
$$\lambda(L, 0) = (-1)^{L/2} \frac{(L-1)!!}{L!!} \quad (L \text{ even})$$

with recurrence for higher M.

### Doubling for $M_T \neq 0$
Each $(M_T, M_U)$ term with $M_T \neq 0$ appears twice (±$M_T$), giving factor of 2.

---

## 7. Transfer Integral (INELDC/GRDSET)

### PRIOR form vertex
$$I_{L_i L_o L_x}^{J_\pi J_\pi'} = \int_0^\infty dr_i \int_0^\infty dr_o \int_{-1}^{1} d(\cos\phi_{ab})\ \chi_{L_i}^{J_\pi}(r_i) \cdot \mathcal{K}(r_i, r_o, \phi_{ab}) \cdot \chi_{L_o}^{J_\pi'*}(r_o)$$

where the kernel is:
$$\mathcal{K} = \phi_T(r_x) \cdot V_{nA}(r_x) \cdot \phi_P(r_p) \cdot \mathcal{A}_{12}(r_i, r_o, \phi_{ab})$$

### POST form vertex (Ptolemy default)
$$\mathcal{K}_{POST} = \phi_T(r_x) \cdot \phi_P(r_p) \cdot V_{np}(r_p) \cdot \mathcal{A}_{12}$$

Note: POST form integrand is large only near $r_p \to 0$ (i.e., $\cos\phi_{ab} \to +1$).

### (U,V) coordinate transformation
$$U = \frac{r_i + r_o}{2}, \quad V = r_i - r_o, \quad \text{Jacobian} = \frac{\partial(r_i,r_o)}{\partial(U,V)} = 1$$

Domain: triangular — $U \in [0, U_{max}]$, $V \in [-2U, +2U]$ (ensures $r_i > 0$, $r_o > 0$).

---

## 8. CUBMAP Adaptive Quadrature (GRDSET)

Ptolemy uses cubic/rational-sinh mapped GL points to concentrate quadrature near the integrand peak.

### CUBMAP formula (rational-sinh, MAPTYP=2)
Given GL nodes $t_k \in [-1,1]$ with weights $w_k$:

$$\tau = \text{arcsinh}(\gamma) = \ln\!\left(\gamma + \sqrt{\gamma^2+1}\right)$$

$$A = -x_{mid}\cdot L, \quad B = L, \quad C = x_{mid}(x_{lo}+x_{hi}) - 2x_{lo}x_{hi}, \quad D = x_{lo}+x_{hi} - 2x_{mid}$$

$$x_k = \frac{-A + (C/\gamma)\sinh(\tau t_k)}{B - (D/\gamma)\sinh(\tau t_k)}, \quad w_k^{new} = w_k \cdot \frac{\tau}{\gamma}\cosh(\tau t_k) \cdot \frac{BC - AD}{[B-(D/\gamma)\sinh(\tau t_k)]^2}$$

Points pile up near $x_{mid}$.

### Default parameters (Ptolemy source.mor line 12300)
| Axis | Map type | $\gamma$ | Midpoint |
|---|---|---|---|
| Sum U=(ri+ro)/2 | MAPTYP=2 (rational-sinh) | 1.0 | SUMMID (adaptive) |
| Dif V=ri-ro | MAPTYP=1 (cubic-sinh) | 5.0 | 0 (symmetric) |
| Phi φ_ab | MAPTYP=2 (rational-sinh) | ~0 (≈linear) | PHI0/2 |

### Adaptive phi cutoff PHI0 per (ri,ro)

For each (ri,ro) pair, Ptolemy does a **2-pass** scan:

**Pass 1:** Scan $\cos\phi$ from +1 downward using test points $x_k = 1 - \frac{2(k-1)^2}{N_{scan}^2}$ (quadratic spacing, dense near $\cos\phi=+1$). Find IEND = last index where $|\text{BSPROD}| > \text{RVRLIM}/(r_i \cdot r_o)$.

**Pass 2:** Set $\phi_0 = \arccos\!\left(1 - \frac{2(\text{IEND}-1)^2}{N_{scan}^2}\right)$, then integrate $\phi \in [0, \phi_0]$ with NPPHI GL points mapped by CUBMAP with midpoint $\phi_0/2$.

$$\text{DPHI}_k = \phi_0 \cdot w_k^{phi} \cdot \sin(\phi_k)$$

This ensures all integration points land where $\mathcal{K}$ is non-negligible.

---

## 9. Transfer S-matrix (SFROMI)

### Assembly
$$S_{SFROMI}(L_i, L_o, L_x, J_\pi, J_\pi') = \text{FACTOR} \cdot \text{ATERM} \cdot \frac{i^{L_i+L_o+2L_x+1}}{\sqrt{2L_i+1}} \cdot I_{L_i L_o L_x}^{J_\pi J_\pi'}$$

### FACTOR
$$\text{FACTOR} = 2\sqrt{\frac{k_a k_b}{E_{cm}^a E_{cm}^b}} = 0.2220 \ \text{MeV}^{-1/2}$$

### ATERM (spectroscopic amplitude)
$$\text{ATERM}(L_x) = \sqrt{\frac{2j_B+1}{2j_A+1}} \cdot \sqrt{2L_x+1} \cdot \mathcal{S}_{AV18} \cdot \mathcal{S}_{target} \cdot W(j_B, l_T, j_A, l_P; j_n, L_x)$$

where:
- $j_B = 1/2$ (neutron in deuteron, 0s1/2) → $2j_B+1 = 2$
- $j_A = 3/2$ (neutron in 34Si, 0d3/2) → $2j_A+1 = 4$
- $W(\ldots)$ = Racah coefficient
- $\mathcal{S}_{AV18} = 0.97069$ = AV18 spectroscopic amplitude for deuteron

### Phase sign
If $(L_i + L_o + 2L_x + 1) \bmod 4 \geq 2$: flip sign of real part.
If $(L_i + L_o + 2L_x + 1) \bmod 2 = 1$: multiply by $i$.

### 9-J structure for J-split
$$S_{SFROMI} = \text{FACTOR} \cdot \mathcal{N}_{9J} \cdot \text{SAV9J} \cdot \text{TEMP} \cdot I$$

where SAV9J and 9J symbols couple the spin quantum numbers of incoming/outgoing channel J-values.

---

## 10. BETCAL: Partial Wave Amplitudes

Sum over incoming partial waves Li:

$$\beta(L_o, L_x, M_x) = \text{FACTOR}_{BET} \cdot \sum_{L_i} (2L_i+1) \cdot C^{L_i\ L_x\ L_o}_{0\ M_x\ M_x} \cdot i^{L_i+L_o+2L_x+1} \cdot e^{2i\sigma_{L_i}^a + 2i\sigma_{L_o}^b} \cdot S_{SFROMI}(L_i, L_o, L_x)$$

where:
$$\text{FACTOR}_{BET} = \frac{1}{2k_a} = 0.3815 \ \text{fm}$$

and $C^{L_i\ L_x\ L_o}_{0\ M_x\ M_x}$ is the Clebsch-Gordan coefficient.

### Legendre normalization for $M_x > 0$
$$\beta(L_o, M_x) \leftarrow \beta(L_o, M_x) \cdot \sqrt{\frac{(L_o - M_x)!}{(L_o + M_x)!}}$$

---

## 11. AMPCAL: Scattering Amplitude

$$F^{M_x}(\theta) = \sum_{L_o} \beta(L_o, L_x, M_x) \cdot P_{L_o}^{M_x}(\cos\theta)$$

where $P_L^M$ = associated Legendre polynomial **with** Condon-Shortley phase (same as `std::assoc_legendre`).

---

## 12. XSECTN: Cross Section

$$\frac{d\sigma}{d\Omega}(\theta) = 10 \cdot \sum_{L_x, M_x} f_{M_x} \cdot |F^{M_x}(\theta)|^2$$

where $f_{M_x} = 1$ for $M_x=0$, $f_{M_x} = 2$ for $M_x > 0$ (sum over $\pm M_x$), and the factor 10 converts fm² to mb.

No AJACOB or JACOB factor appears in the final cross section (absorbed into normalization).

---

## 13. Reference Cross Section (Ptolemy output)

| θ (deg) | dσ/dΩ (mb/sr) |
|---|---|
| 0 | 1.863 |
| 5 | 1.905 |
| 10 | 2.167 |
| 15 | 2.535 |
| 20 | 2.457 |
| 25 | 1.759 |
| 30 | 0.905 |

**Key S-matrix reference values (Ptolemy):**
- Incoming L=0: |S|=0.1496, L=7: 0.6504, L=12: ~1.0
- Outgoing L=0: |S|=0.3716, L=7: 0.9741, L=12: 0.9999
- Transfer (Li=Lo=2, Lx=2): |S_sfromi|=0.01903

---

## 14. Current C++ Status (2026-03-14)

| Component | Status | Notes |
|---|---|---|
| Optical potential | ✅ | WS + SO + Coulomb |
| RCWFN (Coulomb) | ✅ | Steed's method, 6 cases validated |
| WAVELJ (distorted waves) | ✅ | Numerov, S-matrix correct |
| BOUND (bound states) | ✅ | 3 cases validated |
| Coordinate geometry | ✅ | S1,T1,S2,T2 correct |
| A12 angular coupling | ✅ | xlam, 3j, doubling all correct |
| BSPROD clipping | ✅ | phi×V product clipping implemented |
| (U,V) quadrature | ✅ | GL, Jacobian=1, linear interp OK |
| SFROMI / BETCAL | ✅ | ATERM=sqrt(2/4)×..., phases correct |
| AMPCAL / XSECTN | ✅ | PLM with CS phase, factor 10 |
| CUBMAP adaptive phi | 🔄 | **In progress** — needed for POST form |
| POST form | ⏳ | Blocked on CUBMAP |
| 0° cross section | 1.764 mb/sr | 95% of Ptolemy (PRIOR form) |
| Angular shape | ❌ | Falls; should peak at 15° |
