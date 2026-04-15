# Epsilon Algorithm for Partial Wave Summation

## 1. The Problem

In nuclear scattering, the differential cross section involves a partial wave expansion:

$$
f(\theta) = \frac{1}{2ik} \sum_{L=0}^{\infty} (2L+1)\, \beta_L\, P_L(\cos\theta)
$$

where $\beta_L = S_L - 1$ (elastic) or $\beta_L = T_L$ (transfer), and $P_L$ are Legendre polynomials.

In practice we truncate at some $L_{\max}$:

$$
f_{L_{\max}}(\theta) = \frac{1}{2ik} \sum_{L=0}^{L_{\max}} (2L+1)\, \beta_L\, P_L(\cos\theta)
$$

**The convergence problem:** At backward angles ($\theta > 90°$), the Legendre polynomials
$P_L(\cos\theta)$ oscillate with magnitude $\sim 1/\sqrt{L}$ and do not decay. Although $|\beta_L| \to 0$
for $L \gg kR$ (beyond the grazing angular momentum), the Coulomb interaction makes
this decay slow. The partial sums $S_N = \sum_{L=0}^{N}$ oscillate and converge poorly,
requiring very large $L_{\max}$ for stable results at backward angles.

## 2. Sequence Acceleration

Given a sequence of partial sums $\{S_n\}$ that converges to a limit $S$,
a **sequence transformation** maps $\{S_n\}$ to a new sequence $\{T_n\}$
that (ideally) converges faster to the same limit.

Three closely related methods are widely used:

1. **Shanks transformation** (1955)
2. **Wynn's epsilon algorithm** (1956) — an efficient implementation of iterated Shanks
3. **Padé approximants** — the underlying rational approximation theory

### 2.1 Shanks Transformation

The **Shanks transformation** $e_1(S_n)$ is defined as:

$$
e_1(S_n) = \frac{S_{n+1} S_{n-1} - S_n^2}{S_{n+1} - 2S_n + S_{n-1}}
$$

**Derivation:** Assume the error in $S_n$ is approximately geometric:

$$
S_n \approx S + a \cdot r^n
$$

for some constants $a$ and $r$ with $|r| < 1$. Then:

$$
S_{n+1} - S_n = a r^n (r - 1)
$$
$$
S_n - S_{n-1} = a r^{n-1} (r - 1)
$$

Dividing:

$$
\frac{S_{n+1} - S_n}{S_n - S_{n-1}} = r
$$

Solving for $S$ from $S_n = S + a r^n$ and $S_{n+1} = S + a r^{n+1}$:

$$
S = \frac{S_{n+1} S_{n-1} - S_n^2}{S_{n+1} - 2S_n + S_{n-1}} = e_1(S_n)
$$

This is exact when the error is a single geometric term. For more complex
error behavior $S_n \approx S + \sum_j a_j r_j^n$, we iterate: apply Shanks
repeatedly to get $e_2, e_3, \ldots$, each eliminating one more geometric component.

### 2.2 Wynn's Epsilon Algorithm

Computing iterated Shanks transformations directly requires storing many intermediate
sequences and involves numerically unstable determinant ratios. **Wynn (1956)** discovered
an elegant recursive algorithm that computes all iterated Shanks transforms using only
a simple recurrence.

**Definition:** Build a two-dimensional table $\varepsilon_k^{(n)}$ where:
- $n$ = row index (starting partial sum)
- $k$ = column index (level of acceleration)

**Initialization:**

$$
\varepsilon_{-1}^{(n)} = 0, \qquad \varepsilon_0^{(n)} = S_n
$$

**Recurrence (the rhombus rule):**

$$
\varepsilon_{k+1}^{(n)} = \varepsilon_{k-1}^{(n+1)} + \frac{1}{\varepsilon_k^{(n+1)} - \varepsilon_k^{(n)}}
$$

The table looks like:

```
n\k    -1        0           1              2
 0      0       S_0
                     1/(S_1-S_0)
 1      0       S_1                     e_1(S_1)
                     1/(S_2-S_1)
 2      0       S_2                     e_1(S_2)
                     1/(S_3-S_2)
 3      0       S_3                     e_1(S_3)
```

**Key property:** The even columns give the Shanks transforms:

$$
\varepsilon_{2p}^{(n)} = e_p(S_n)
$$

The odd columns are intermediate reciprocal differences (not directly useful).

**Why it works:** Each application of the rhombus rule effectively performs one level of
Shanks transformation. The recurrence avoids the explicit determinant computation,
making it numerically stable and efficient — O(N²) operations for N partial sums,
with only O(N) storage if computed column-by-column.


### 2.3 Connection to Padé Approximants

A **Padé approximant** $[M/N]$ to a power series $\sum_{k=0}^{\infty} c_k x^k$ is a
rational function $P_M(x)/Q_N(x)$ (polynomial of degree $M$ over polynomial of degree $N$)
that matches the first $M + N + 1$ coefficients of the series.

The connection to partial wave sums is through the **formal series:**

$$
S(x) = \sum_{n=0}^{\infty} a_n x^n
$$

evaluated at $x = 1$. The partial sums $S_n = \sum_{k=0}^{n} a_k$ are the
values of the truncated polynomial at $x = 1$.

**Theorem (Wynn, 1956):** The even-column entries of the epsilon table satisfy:

$$
\varepsilon_{2p}^{(n)} = [n + p \,/\, p] \text{ Padé approximant of } S(x) \text{ at } x=1
$$

In other words, the epsilon algorithm implicitly constructs a sequence of diagonal and
near-diagonal Padé approximants of increasing order, evaluated at the point of interest.

**Why Padé helps for partial waves:** A polynomial (truncated Legendre sum) has no
poles and cannot represent the rapidly varying Coulomb amplitude efficiently. A
rational function (Padé) can place poles to mimic the branch-cut structure of the
Coulomb amplitude, giving much faster convergence with fewer terms.

### 2.4 Convergence Properties

For a sequence with error of the form:

$$
S_n - S = \sum_{j=1}^{p} a_j r_j^n + O(r_{p+1}^n)
$$

the $p$-th Shanks transform $e_p$ (equivalently $\varepsilon_{2p}$) eliminates
all $p$ geometric components exactly, leaving an error of order $r_{p+1}^n$.

In the partial wave context, each geometric component corresponds to a
diffraction/interference pattern. The epsilon algorithm can "see through"
the oscillations and extract the converged limit, even when the raw partial
sums are still oscillating wildly.

**Limitations:**
- Breaks down if $\varepsilon_k^{(n+1)} \approx \varepsilon_k^{(n)}$ (division by near-zero)
- Not guaranteed to converge for all series (but works well for the oscillatory series typical in scattering)
- Complex arithmetic required for complex amplitudes (Ptolemy uses complex epsilon)


## 3. Implementation in Fortran Ptolemy

In `source.f`, subroutine `AMPCAL` (line ~220) computes the scattering amplitude
at each angle. The key section (lines ~470-500):

```fortran
C     IF NECESSARY APPLY EPSILON ALGORITHM
C
 150     FERROR(KOFFS) = ( DABS(CONTR(LOMX+1))+DABS(CONTI(LOMX+1)) )
     1      / ( DABS(FT(1))+DABS(FT(2)) + SMLNUM )
         IF ( LEBACK .LE. 0 )  GO TO 200
         N = MIN0( (LOMX-LOMNMX)/LSKP, LEBACK-1 )
         IF ( N .LE. 5 ) GO TO 200
         ...
         CALL EPSLON ( FEPSLO, NN, FT, 1.D-5, FERROR(KOFFS), IER )
```

The algorithm:
1. Accumulates partial sums as each $L$ term is added
2. Stores the last `LEBACK` partial sums in array `FEPSLO`
3. Calls `EPSLON` to apply Wynn's algorithm to these partial sums
4. Replaces the raw sum `FT` with the epsilon-accelerated result
5. Reports the estimated convergence error in `FERROR`

The `EPSLON` subroutine (line ~15261) implements the complex version of
Wynn's algorithm with convergence checking at tolerance $10^{-5}$.

## 4. Implications for C++ Ptolemy++

The current C++ elastic solver (`elastic.cpp`) uses a **direct Legendre sum**
without epsilon acceleration. This means:

- **Forward angles** ($\theta < 40°$): Good convergence, no acceleration needed.
  The partial wave sum is dominated by high-$L$ terms that converge rapidly.
- **Backward angles** ($\theta > 90°$): Poor convergence, large truncation errors.
  The epsilon algorithm is essential for reliable results without extremely large $L_{\max}$.

**Recommended implementation:** Add Wynn's epsilon algorithm as a post-processing
step on the partial wave sum in `DCSUnpolarized()`. The algorithm is simple (~50 lines)
and requires only O($L_{\max}$) additional storage.

## 5. Worked Example

Consider a geometric sequence $S_n = 1 + r + r^2 + \cdots + r^n$ converging to $1/(1-r)$.

For $r = 0.9$: $S_{10} = 6.8619$, true limit = $10.0$, error = 31%.

Applying one Shanks transform:

$$
e_1(S_{10}) = \frac{S_{11} \cdot S_9 - S_{10}^2}{S_{11} - 2 S_{10} + S_9}
= \frac{7.1757 \times 6.5132 - 6.8619^2}{7.1757 - 2 \times 6.8619 + 6.5132}
= \frac{46.733 - 47.086}{-0.035} = 10.08
$$

One Shanks iteration reduced the error from 31% to 0.8% — because the error
was almost exactly geometric ($a \cdot r^n$ with $r = 0.9$).

For the partial wave sum, the error structure is more complex (multiple geometric
components from different diffraction patterns), but iterated Shanks (= epsilon
algorithm) still provides dramatic acceleration.

## References

1. D. Shanks, "Non-linear transformations of divergent and slowly convergent sequences,"
   *J. Math. Phys.* **34**, 1-42 (1955).
2. P. Wynn, "On a device for computing the $e_m(S_n)$ transformation,"
   *Math. Tables Aids Comput.* **10**, 91-96 (1956).
3. C. Brezinski and M. Redivo-Zaglia, "Extrapolation Methods: Theory and Practice,"
   North-Holland (1991).
4. M.H. Macfarlane and S.C. Pieper, "Ptolemy: A Program for Heavy-Ion Direct-Reaction
   Calculations," ANL-76-11 (1978).

