# Bound State Calculation Algorithm

This document describes the numerical method used to calculate the nuclear bound state wave functions and potentials in the C++ DWBA implementation. The algorithm replicates the robust "Matching Method" used in the original Ptolemy Fortran code.

## Problem Statement
We solve the radial Schrödinger equation for a bound state:
$$ \left[ -\frac{\hbar^2}{2\mu} \frac{d^2}{dr^2} + V_{eff}(r) \right] u(r) = E u(r) $$
where $E < 0$ is the binding energy. The effective potential $V_{eff}(r)$ includes the nuclear potential (Woods-Saxon), spin-orbit coupling, Coulomb potential, and the centrifugal barrier.

The goal is to find the potential depth $V_0$ (or other parameters) such that a solution exists with the correct binding energy $E$ and the specified number of nodes $n$ (excluding the origin).

## Algorithm: Inward-Outward Matching Method

Unlike simple shooting methods that integrate from the origin to infinity (which are unstable for bound states due to the exponentially growing solution in the forbidden region), this method integrates from both sides and matches them at a suitable radius.

### 1. Grid and Matching Radius
- **Grid**: Uniform mesh with step size $h$ (typically 0.1 fm).
- **Matching Radius ($R_{match}$)**: Chosen to be slightly outside the nuclear radius, typically $R_{match} \approx R_0 A^{1/3} + 2.0$ fm. This ensures the matching occurs in a region where the potential is changing but not yet purely asymptotic.

### 2. Outward Integration ($0 \to R_{match}$)
- **Boundary Condition**: $u(0) = 0$, $u(h) \sim h^{l+1}$.
- **Method**: Numerov integration.
- **Node Counting**: The number of nodes is counted during integration.
- **Result**: $u_{out}(r)$ defined for $r \in [0, R_{match}]$.

### 3. Inward Integration ($R_{max} \to R_{match}$)
- **Boundary Condition**: At large $R_{max}$, the solution decays exponentially:
  $$ u(R_{max}) \sim \exp(-\kappa R_{max}) $$
  where $\kappa = \sqrt{2\mu |E|} / \hbar$.
- **Method**: Numerov integration (backward).
- **Result**: $u_{in}(r)$ defined for $r \in [R_{match}, R_{max}]$.

### 4. Matching Condition
At $R_{match}$, the two solutions must match in value and slope (logarithmic derivative).
- **Scaling**: Scale $u_{in}$ such that $u_{in}(R_{match}) = u_{out}(R_{match})$.
- **Logarithmic Derivative**:
  $$ \text{Diff} = \left( \frac{u'}{u} \right)_{out} - \left( \frac{u'}{u} \right)_{in} $$
  The derivatives are calculated using the Numerov relation or 3-point formulas consistent with the grid.

### 5. Iteration Strategy
We iterate on the potential depth $V$ to drive $\text{Diff} \to 0$.

1.  **Node Check**:
    - If `nodes > n`: Potential is too deep. Decrease $V$.
    - If `nodes < n`: Potential is too shallow. Increase $V$.
    - Step size is reduced as we approach the correct node count.

2.  **Newton-Raphson / Gradient Descent**:
    - Once the node count is correct, we use the matching difference `Diff` to update $V$.
    - Update rule: $V_{new} = V_{old} + \text{factor} \times \text{Diff}$.
    - **Adaptive Factor**: The `factor` starts at 2.0 (empirically matched to Fortran behavior). If oscillation is detected (sign of `Diff` flips), the factor is reduced by half.

## Verification
This algorithm has been verified against the original Fortran `BOUND` subroutine.
- **Test Case 1 (Target)**: Converges to $V \approx 47.07$ MeV (Matches Fortran).
- **Test Case 2 (Projectile)**: Converges to $V \approx 62.91$ MeV (Matches Fortran).

The consistency between C++ and Fortran implementations confirms the correctness of the solver logic.
