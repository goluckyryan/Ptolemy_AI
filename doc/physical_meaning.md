# Scattering Theory — Physical Meaning
*Session: 2026-04-19, Ryan × Dudu*

---

## 1. Elastic Scattering — the standard picture

The full wavefunction has the asymptotic form:

$$\Psi \xrightarrow{r\to\infty} e^{ikz} + f(\theta)\frac{e^{ikr}}{r}$$

Steps:
1. **Solve the SE** for each partial wave (L, J) via Numerov outward integration, with optical potential V_C + U_N
2. **Extract S_LJ** by matching to Coulomb functions F_L, G_L at large r (Wronskian or 2-point matching)
3. **Build the nuclear amplitude** as a sum over L:

$$f_N(\theta) = \frac{1}{2ik}\sum_L (2L+1)\,e^{2i\sigma_L}(S_L - 1)\,P_L(\cos\theta)$$

4. **Add Rutherford amplitude** f_C(θ) (closed-form, pure Coulomb)
5. **DCS** = |f_C + f_N|² / (2S+1), averaged over initial spin projections

Key: each term carries **e^{2iσ_L}** (Coulomb phase shift σ_L = arg Γ(L+1+iη)) — without these, the interference between f_C and f_N is wrong.

For spin > 0 (proton, deuteron): the amplitude becomes a **(2S+1)×(2S+1) matrix** in spin space, with spin-no-flip (v=v₀, uses P_L) and spin-flip (v≠v₀, uses P_L^{|Δm|}) terms. The Coulomb amplitude only contributes to spin-no-flip.

---

## 2. Why Partial Waves are a Shortcut

The T-matrix (integral form) for elastic:

$$f_N(\theta) = -\frac{\mu}{2\pi\hbar^2}\langle \chi^{(-)}_{k'} | U_N | \psi^{(+)}_k \rangle$$

The partial wave sum is derived from this by:
1. Expanding χ^{(±)} in partial waves (each carrying e^{iσ_L})
2. Applying **Green's theorem** — the volume integral collapses to a **surface term** at r→∞
3. The surface term evaluates via the Wronskian to exactly (S_L − 1) × kinematic factor

So the partial wave sum B is derived directly from the T-matrix A — they are the same formula in different languages.

**The shortcut:** Partial waves exploit spherical symmetry to replace one 3D problem with N independent 1D problems. The "integration" collapses to reading S_L from the asymptotics of u_L(r). No explicit integral needed.

This only works because:
- The potential is **central** (L is a good quantum number)
- The same potential that scatters also defines the asymptotic channel → Green's theorem connects interior to surface

---

## 3. The Physics of the T-matrix

### Pure Coulomb case

The Coulomb wave χ is already the **exact eigenstate** of H_C. The Rutherford amplitude comes entirely from reading the asymptotic tail of χ — the outgoing e^{ikr}/r part with coefficient f_C(θ). No integral needed. T_C = 0 for k ≠ k' (eigenstate orthogonality).

### With nuclear potential

ψ ≠ χ in the interior (nuclear range, r < few fm). The T-matrix:

$$T_N = \langle \chi^{(-)}_{k'} | U_N | \psi^{(+)} \rangle$$

**Physical meaning:**
- U_N acts only in the interior where it's nonzero
- χ^{(-)} acts as the **detector** — the natural basis for "what comes out" in the Coulomb field
- U_N is the **cause** of the deviation from pure Coulomb
- The integral picks up only where U_N ≠ 0 (the nuclear surface region)

The effect of U_N is already encoded in ψ^{(+)} — the T-matrix is just a way to extract it. For elastic, T = ⟨k'|U|ψ⟩ is mathematically equivalent to "read f(θ) from the tail of ψ directly." U is redundant for elastic — the information is already in the asymptotics.

---

## 4. How the Interior Integral Equals the Asymptotic Tail (Green's Theorem)

This is the key bridge between the T-matrix picture and the partial wave picture.

Both χ_L and ψ_L satisfy Schrödinger equations:
- χ_L: (T + V_C − E) χ_L = 0
- ψ_L: (T + V_C + U_N − E) ψ_L = 0

Subtract one from the other:

$$\chi_L^*(T - E)\psi_L - \psi_L(T - E)\chi_L^* = \chi_L^* U_N \psi_L$$

The left side is a **perfect divergence** (kinetic energy is a second derivative). Integrate over all space:

$$\underbrace{\int_0^\infty \chi_L^* U_N \psi_L\, dr}_{\text{interior T-matrix integral}} = \underbrace{\Big[\chi_L^* \psi_L' - \psi_L \chi_L^{*\prime}\Big]_0^\infty}_{\text{surface (Wronskian) at } r\to\infty}$$

**The interior integral has completely disappeared — replaced by a boundary term at infinity.**

Inserting the asymptotics at large r (χ_L → F_L, ψ_L → (H_L^− − S_L H_L^+)/2i) and evaluating the Wronskians:

$$W[F_L, H_L^+] = -\frac{1}{k}, \qquad W[F_L, H_L^-] = +\frac{1}{k}$$

gives:

$$\int_0^\infty \chi_L^* U_N \psi_L\, dr \;\propto\; \frac{S_L - 1}{k}$$

**The interior integral = (S_L − 1)/k. This is exactly what the partial wave sum uses.**

### Physical analogy: Gauss's Law

Think of it like electrostatics:
- **Source** = U_N in the interior (the nuclear potential region, r < few fm)
- **Flux** = the Wronskian at the surface (the asymptotic tail of ψ)

The nuclear potential modifies ψ only in the interior. But that modification **propagates outward** and shows up at infinity as the deviation of S_L from pure Coulomb (i.e. S_L ≠ e^{2iσ_L}). Green's theorem is the statement that **what happens inside is perfectly encoded in what comes out at the surface** — because ψ is continuous and smooth everywhere in between.

The interior integral and the asymptotic tail are **not two different things** — they are the same information read from different places. Green's theorem is the pipeline connecting them.

---

## 5. Why DWBA Requires a Genuine Integral

For transfer reactions (e.g. d+A → p+B):

- Initial channel: d+A, distorted wave χ^{(+)}_i generated by U_optical(d+A)
- Final channel: p+B, distorted wave χ^{(-)}_o generated by U_optical(p+B)
- Transfer interaction: V_np (neutron-proton force inside the deuteron)

**V_np is genuinely new** — it was never used to build either distorted wave. It's the physics that **causes the transition** between channels. Unlike the elastic case, you cannot read it off from the asymptotics of χ_i or χ_o, and Green's theorem cannot be applied (V_np is not the potential of either channel).

The DWBA T-matrix:

$$T_{fi} = \langle \chi^{(-)}_{k_o}\phi_p\phi_B \;|\; V_{np}\cdot\phi_P(r_p) \;|\; \phi_A\phi_T(r_T)\chi^{(+)}_{k_i}\rangle$$

**Physical picture:** Two distorted waves meet in the nuclear interior. V_np causes the transfer — it rips the neutron out of the deuteron (φ_P = projectile bound state) and deposits it into the target (φ_T = target bound state). The overlap integral measures how effectively that kick connects the incoming to the outgoing channel.

The integral is unavoidable because:
1. V_np is **not encoded** in either χ — Green's theorem does not apply
2. Initial and final states are **different particles** — you must project explicitly onto the outgoing channel
3. DWBA factorizes ψ^{(+)} ≈ χ_i · φ_d · φ_A, which forces V_np back in as an explicit integral

**The U in the DWBA T-matrix is the price of the DWBA approximation.** If you knew the exact many-body ψ^{(+)}, you could read off the transfer amplitude from the asymptotics of the p+B channel — but you don't, so you integrate.

---

## 6. Comparison Table

| | Elastic | DWBA Transfer |
|---|---|---|
| In/out channels | same (d+A → d+A) | different (d+A → p+B) |
| Sum variable | single L | double sum (Li, Lo) |
| "S-matrix" | S_L = e^{2i(σ+δ)}, scalar | T_{Li,Lo} = radial integral, complex matrix |
| Angular function | P_L(cosθ) | CG × P_{Lx}^m(cosθ) |
| T-matrix integral | collapses to surface term = asymptotics | genuine volume integral, no shortcut |
| Green's theorem applies? | ✅ yes — U_N built ψ | ❌ no — V_np not in either χ |
| New interaction | U_N already in ψ (redundant in T) | V_np genuinely new (essential in T) |
| Code | Numerov → match → S_L → sum | Numerov × 2 → explicit radial integral (INRDIN/GRDSET) |

---

## 7. TODO: Polarization Observables

For spin-S projectile, the full amplitude is a (2S+1)×(2S+1) matrix:

$$f_{m_s', m_s}(\theta) = f_C(\theta)\,\delta_{m_s'm_s} + f_N^{m_s'm_s}(\theta)$$

Observables to implement:
- **Analyzing power** A_y: from Tr[F σ_y F†] / Tr[F F†]
- **Tensor analyzing powers** iT11, T20, T21, T22 (deuteron S=1)
- Requires storing the full amplitude matrix F at each angle

Currently `elastic.cpp` computes only DCSUnpolarized. Polarization is on the TODO list.

---

## 8. Zero-Range Approximation (ZR)

*To be discussed in a future session.*

The ZR approximation replaces the 6D DWBA integral with a 1D integral:

$$V_{np}(r_p)\,\phi_P(r_p) \approx D_0\,\delta^3(r_p)$$

with D₀ ≈ −120.1 MeV·fm^{3/2}. This collapses r_p → 0, leaving only the target coordinate r_a = r_b (up to a kinematic scale factor zr_scale = A/(A+1)). The 6D integral becomes 1D.

---

## 9. Why f(θ) Has No r Dependence

The scattered wave e^{ikr}/r is the **asymptotic form** — valid only at large r, far outside the potential. In that region the potential is zero, no force is acting, and the particle is free.

The free Schrödinger equation at large r:

$$\nabla^2\psi + k^2\psi = 0$$

The general outgoing spherical wave solution is:

$$\psi_{\text{scattered}} = \frac{g(\theta, r)}{r}\,e^{ikr}$$

Substituting into the free SE, the radial equation for g at large r gives:

$$g(r) \to \text{constant as } r\to\infty$$

That constant can depend on θ (set by what happened in the scattering region) but **not on r** — because the particle is free and simply propagating outward. That constant is f(θ).

### The 1/r factor — flux conservation

The 1/r is not mysterious — it's **flux conservation**. The same total probability must pass through every sphere of radius r:

$$|\psi|^2 \times 4\pi r^2 = \text{const} \implies |\psi| \propto \frac{1}{r}$$

So amplitude must fall as 1/r as the wave spreads outward.

### Physical picture

Think of ripples on water. A stone dropped in the center creates ripples spreading outward. Far from the stone:
- The **shape** (angular pattern) of each ripple is fixed — determined by what happened at the center
- As the ripple travels outward it weakens (1/r in amplitude) but the angular pattern f(θ) does not change

f(θ) is **frozen** once the wave exits the interaction region. The r-dependence is just trivial outward spreading.

### Summary
- **e^{ikr}/r** = outgoing spherical wave, flux conservation
- **f(θ)** = angular pattern set by the scattering potential, frozen once the wave exits the potential region, independent of r

---

*Notes by Dudu 🐱 — 2026-04-19*
