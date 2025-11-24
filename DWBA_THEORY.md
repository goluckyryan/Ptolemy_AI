# DWBA Theory and Ptolemy Implementation

This document outlines the theoretical background of the Distorted Wave Born Approximation (DWBA) and details how the Ptolemy code (and our C++ translation) implements the **Finite-Range** version of this theory. It contrasts this with the simpler **Zero-Range** approximation often found in introductory handbooks.

## 1. Elastic Scattering (The Foundation)

Before calculating transfer reactions, we must describe the elastic scattering of the projectile and ejectile in the entrance and exit channels. This provides the "Distorted Waves" ($\chi$).

### Physics
The relative motion of two nuclei is governed by the Schrödinger equation with an **Optical Potential** $U(r)$:
$$ \left[ -\frac{\hbar^2}{2\mu} \nabla^2 + U(r) - E \right] \chi(\vec{r}) = 0 $$

$U(r)$ contains:
-   **Real Central**: Nuclear attraction (Woods-Saxon).
-   **Imaginary**: Absorption (Volume/Surface).
-   **Spin-Orbit**: Coupling between orbital angular momentum $L$ and spin $S$.
-   **Coulomb**: Repulsion between protons.

### Partial Wave Expansion
The wave function $\chi(\vec{r})$ is expanded in partial waves:
$$ \chi(\vec{r}) = \frac{4\pi}{k} \sum_{L,M} i^L \frac{u_L(r)}{r} Y_{LM}(\hat{r}) Y_{LM}^*(\hat{k}) $$

The radial function $u_L(r)$ satisfies:
$$ \frac{d^2 u_L}{dr^2} + \left[ k^2 - \frac{2\mu}{\hbar^2} U(r) - \frac{L(L+1)}{r^2} \right] u_L(r) = 0 $$

### Implementation ([WavElj](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#394-554))
-   **Algorithm**: Numerov method (4th order) or Runge-Kutta to integrate $u_L(r)$ from the origin ($r=0$) to a matching radius ($R_{match}$).
-   **Matching**: At $R_{match}$ (where nuclear potential is negligible), the internal solution is matched to the asymptotic Coulomb functions ($F_L, G_L$) to determine the **Phase Shift** $\delta_L$ and **S-Matrix** $S_L$.
    $$ u_L(r) \xrightarrow{r \to \infty} \frac{i}{2} [ H_L^-(kr) - S_L H_L^+(kr) ] $$
-   **Code**: `DWBA::WavElj` performs this integration and calculates $S_L$.

---

## 2. Transfer Reaction $A(a,b)B$

We consider a reaction where a particle $x$ is transferred from projectile $a$ to target $A$:
$$ a (= b + x) + A \to b + B (= A + x) $$

### The Transition Amplitude ($T_{ba}$)
In DWBA, the transition amplitude is:
$$ T_{ba} = \int d^3r_b \int d^3r_a \chi_b^{(-)*}(\vec{r}_b) \langle \psi_B \psi_b | V | \psi_A \psi_a \rangle \chi_a^{(+)}(\vec{r}_a) $$

-   $\chi_a, \chi_b$: Distorted waves for entrance ($a+A$) and exit ($b+B$) channels.
-   $\psi$: Internal wave functions.
-   $V$: Interaction potential (usually $V_{bx}$, the potential binding $x$ to $b$ in the prior form).

Separating the internal states (Spectroscopic Amplitudes) from the radial motion:
$$ T_{ba} \propto \int d^3r_b \int d^3r_a \chi_b^*(\vec{r}_b) \phi_{Bx}^*(\vec{r}_x) V(\vec{r}_{bx}) \phi_{ax}(\vec{r}_{bx}) \chi_a(\vec{r}_a) $$

-   $\phi_{Bx}$: Bound state wave function of $x$ in $B$ (Target Bound State).
-   $\phi_{ax}$: Bound state wave function of $x$ in $a$ (Projectile Bound State).
-   $\vec{r}_a, \vec{r}_b$: Channel coordinates.
-   $\vec{r}_x, \vec{r}_{bx}$: Bound state coordinates (depend on $\vec{r}_a, \vec{r}_b$).

---

## 3. Zero-Range Approximation (Handbook Approach)

The **Zero-Range (ZR)** approximation assumes the interaction $V(\vec{r}_{bx})$ and the projectile bound state $\phi_{ax}(\vec{r}_{bx})$ are non-zero only when $b$ and $x$ are at the same point (i.e., $a$ is a point particle).

$$ V(\vec{r}_{bx}) \phi_{ax}(\vec{r}_{bx}) \approx D_0 \delta(\vec{r}_{bx}) $$

This implies $\vec{r}_b \propto \vec{r}_a$. The 6-dimensional integral collapses to a **3-dimensional integral**:
$$ T_{ZR} \propto \int d^3r \chi_b^*(\beta \vec{r}) \phi_{Bx}^*(\vec{r}) \chi_a(\vec{r}) $$

This is computationally very fast but inaccurate for:
-   Heavy ion reactions.
-   High energy reactions.
-   States with non-zero orbital angular momentum transfer where finite size matters.

---

## 4. Finite-Range DWBA (Ptolemy Approach)

Ptolemy performs a **Finite-Range (FR)** calculation, treating the full spatial extent of the wave functions and interaction.

### Coordinate Transformation
We must integrate over two independent vectors, e.g., $\vec{r}_a$ and $\vec{r}_b$. The bound state coordinates are linear combinations:
$$ \vec{r}_x = S_1 \vec{r}_a + T_1 \vec{r}_b $$
$$ \vec{r}_{bx} = S_2 \vec{r}_a + T_2 \vec{r}_b $$

### Partial Wave Expansion
To evaluate the 6D integral, we expand all wave functions in partial waves (Spherical Harmonics).
$$ T \sim \sum_{L_a, L_b} \int r_a^2 dr_a \int r_b^2 dr_b \chi_{L_b}(r_b) K_{L_a L_b}(r_a, r_b) \chi_{L_a}(r_a) $$

The Kernel $K$ involves an angular integration over the angle $\theta$ between $\vec{r}_a$ and $\vec{r}_b$:
$$ K(r_a, r_b) = \int_{-1}^{1} d(\cos \theta) \phi_{Bx}^*(r_x) V(r_{bx}) \phi_{ax}(r_{bx}) P_L(\cos \theta) $$

where $r_x$ and $r_{bx}$ depend on $r_a, r_b, \theta$ via the law of cosines:
$$ r_x^2 = s_1^2 r_a^2 + t_1^2 r_b^2 + 2 s_1 t_1 r_a r_b \cos \theta $$

### Ptolemy Algorithm ([InelDc](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#568-792))
1.  **Grids**: Define grids for $r_a$ (Incoming), $r_b$ (Outgoing), and $\theta$.
2.  **Bound States**: Calculate $\phi_{Bx}$ and $\phi_{ax}$ on these grids (interpolating if necessary).
    -   [CalculateBoundState](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#896-1027): Solves for $\phi$ using the Matching Method.
3.  **Distorted Waves**: Calculate $\chi_{L_a}(r_a)$ and $\chi_{L_b}(r_b)$ for all relevant partial waves.
    -   [WavElj](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#394-554): Solves the optical model equation.
4.  **Integration**:
    -   Loop over $L_a$ (Incoming Partial Wave).
    -   Loop over $L_b$ (Outgoing Partial Wave) allowed by selection rules.
    -   Loop over $r_a, r_b$:
        -   Loop over $\theta$: Calculate form factor $F = \phi_{Bx} V \phi_{ax}$.
        -   Sum angular contribution.
    -   Sum radial contribution weighted by $\chi_a \chi_b$.
5.  **S-Matrix**: The result is the Transfer S-Matrix element $S_{L_a, L_b}$.
6.  **Cross Section**: Sum $|S|^2$ with geometric factors to get $d\sigma/d\Omega$.

### Key Differences from Zero-Range
-   **Complexity**: Double radial integral + Angular integral vs. Single radial integral.
-   **Accuracy**: Correctly handles recoil and finite size of the projectile.
-   **Computational Cost**: Significantly higher ($N_r^2 \times N_\theta$ vs $N_r$).

## 5. Code Mapping

| Theory Component | Ptolemy/C++ Function | Description |
| :--- | :--- | :--- |
| **Bound State $\phi$** | [CalculateBoundState](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#896-1027) | Solves radial Schrödinger eq for bound state (negative energy). Uses Matching Method. |
| **Distorted Wave $\chi$** | [WavElj](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#394-554) | Solves radial Schrödinger eq for scattering state (positive energy). Uses Numerov. |
| **Elastic S-Matrix** | [WavElj](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#394-554) | Extracted from asymptotic matching of $\chi$. |
| **Integration Grid** | [GrdSet](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#555-570) | Sets up $\theta$ grid (Gauss-Legendre). $r$ grids from [WavSet](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#343-383). |
| **Transfer Kernel** | [InelDc](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#568-792) (Inner Loop) | Calculates $r_x, r_{bx}$ and form factor $\phi_T V \phi_P$. |
| **Radial Integral** | [InelDc](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#568-792) (Outer Loops) | Integrates $\chi_b^* \times \text{Kernel} \times \chi_a$. |
| **Cross Section** | [XSectn](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#823-925) | Sums partial wave amplitudes to get differential cross section. |

## 6. Current Status of C++ Implementation
-   **Bound States**: Fully implemented and verified (Matching Method). Matches Ptolemy results.
-   **Distorted Waves**: Implemented (Numerov). Elastic S-Matrix calculation needs verification.
-   **Integration**: Finite-Range integration logic implemented in [InelDc](file:///home/ryan/ptolemy_2019/Cpp_AI/src/dwba/dwba.cpp#568-792).
-   **Cross Section**: Formula implemented, but currently yielding zero (likely normalization/prefactor issue).
