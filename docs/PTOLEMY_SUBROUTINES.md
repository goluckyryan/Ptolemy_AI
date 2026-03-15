# Ptolemy Subroutine Map — DWBA Calculation Flow

> **Reaction context:** 33Si(d,p)34Si at E_lab=20 MeV, L=2 transfer  
> **Target:** Match cross section at 0° = 1.863 mb/sr  
> **Source:** `src/source.mor`, `src/rcwfn.f`, `src/fortlib.mor`

---

## 🗺️ High-Level Flow

```
┌─────────────────────────────────────────────────────────────┐
│                        CONTRL                               │
│            (Main dispatcher / sequencer)                    │
│                    [line 8109]                              │
└──────────┬──────────────────────────────────────────────────┘
           │
     ┌─────▼──────┐
     │  DATAIN    │  Parse input file (potentials, kinematics,
     │ [line 11643]│  bound state params, angles)
     └─────┬──────┘
           │
     ┌─────▼──────┐
     │   BOUND    │  Compute bound state wavefunctions φ(r)
     │ [line 3642] │  for target (34Si) and projectile (deuteron)
     └─────┬──────┘
           │
     ┌─────▼──────┐
     │  WAVSET    │  Set up distorted wave computation
     │ [line 31940]│  (allocate memory, set grid parameters)
     └─────┬──────┘
           │
     ┌─────▼──────┐
     │  GETSCT    │  Compute optical model S-matrices
     │ [line 15538]│  (calls WAVELJ for each L)
     └─────┬──────┘
           │
     ┌─────▼──────┐
     │  GRDSET    │  Set up integration grid (Gauss points)
     │ [line 15710]│  for the radial DWBA integral
     └─────┬──────┘
           │
     ┌─────▼──────┐
     │  INELDC    │  Main DWBA radial integral
     │ [line 17454]│  Computes T-matrix elements I(Li,Lo,Lx)
     └─────┬──────┘
           │  (calls SFROMI inside loop)
     ┌─────▼──────┐
     │  XSECTN    │  Orchestrate cross section output
     │ [line 32743]│  (calls BETCAL then AMPCAL then ANAPOW)
     └─────┬──────┘
           │
     ┌─────▼──────┐
     │  BETCAL    │  Compute angle-independent beta amplitudes
     │ [line 3358] │  β(Lo) from S_sfromi elements
     └─────┬──────┘
           │
     ┌─────▼──────┐
     │  AMPCAL    │  Compute angular distribution F(θ)
     │  [line 220] │  via Legendre polynomial sum
     └─────┬──────┘
           │
     ┌─────▼──────┐
     │  ANAPOW    │  Compute and print d²σ/dΩ vs angle
     │  [line 578] │  (final cross section in mb/sr)
     └────────────┘
```

---

## 📦 Subroutine Reference

### 🔷 CONTRL — Main Dispatcher
**File:** `source.mor` line 8109  
**Role:** Sequences the entire Ptolemy calculation. Acts as a state machine — calls subroutines in order based on `IGOTO` parameter.

| IGOTO | Action |
|-------|--------|
| 1 | Call BOUND (after projectile/target input) |
| 5 | Call WAVSET (after incoming/outgoing setup) |
| 7 | Call INELDC (after RADIALIN grid setup) |
| 9 | Call XSECTN (after CROSSSEC) |

---

### 🔷 DATAIN — Input Parser
**File:** `source.mor` line 11643  
**Role:** Reads the Ptolemy input file. Parses:
- Optical model potentials (Woods-Saxon parameters V, r₀, a for each channel)
- Kinematics (Elab, Q-value, masses)
- Bound state parameters (n, l, j, binding energy)
- Angular range and step

**Outputs:** All parameters stored in COMMON blocks for use by all subroutines.

---

### 🔷 BOUND — Bound State Wavefunction
**File:** `source.mor` line 3642  
**Role:** Computes the radial bound state wavefunction for the transferred nucleon.

**Physics:**
- Solves the Schrödinger equation with a Woods-Saxon + spin-orbit potential
- Uses Numerov integration from r=0 outward
- Performs depth search: adjusts V until binding energy matches
- **Stores:** φ(r) = u(r)/r (divided by r), where ∫u²dr = 1

**Inputs:**
```
n, l, j      - quantum numbers
E_binding    - binding energy (MeV)
V, r₀, a     - WS central potential
Vso, Rso, aso - spin-orbit potential
Rc           - Coulomb radius
A            - mass number
```

**Key detail:** Normalization is ∫φ²r²dr = 1 ↔ ∫u²dr = 1

**Called for:** target bound state (n in ³⁴Si, 0d3/2) AND projectile bound state (n in deuteron, 0s1/2)

---

### 🔷 WAVSET — Distorted Wave Setup
**File:** `source.mor` line 31940  
**Role:** Allocates memory and initializes the distorted wave calculation. Sets grid parameters (step size H, number of steps NSTEP) for WAVELJ.

**Calls:** WAVELJ (via WFGET) for each partial wave L

---

### 🔷 WAVELJ — Distorted Wave (Numerov Integration)
**File:** `source.mor` line 30428  
**Role:** Computes the distorted wave u_L(r) = r·χ_L(r) for a single partial wave L using the modified Numerov algorithm.

**Physics:**
```
[-d²/dr² + L(L+1)/r² + U(r)] u_L(r) = k² u_L(r)

where U(r) = 2μ/ℏ² × [V_WS(r) + V_SO(r) + V_Coulomb(r)]
```

**Algorithm:**
1. Start at r=0 with u_L(0) = 0, u_L(h) = h^(L+1)
2. Integrate outward using W-form modified Numerov
3. Match to Coulomb functions F_L(kr) and G_L(kr) at large r
4. Extract S-matrix: S_L = e^(2iδ_L)

**Key equation (from source.mor):**
```fortran
DLSQ = DL*(DL+1)      ! L*(L+1) centrifugal term
WAVR(I+4) = VREAL(I+1) + DLSQ*VCENT(I+1)   ! effective potential
! VCENT = 1/r²
```

**Output stored as:** u_L(r) = r·χ_L(r)  ← **NOT χ_L(r) directly!**

**Reference S-matrix values (Ptolemy, this reaction):**
| L | |S| incoming (d+³³Si) | |S| outgoing (p+³⁴Si) |
|---|---|---|
| 0 | 0.1496 | 0.3716 |
| 1 | 0.1710 | 0.4582 |
| 2 | 0.1283 | 0.2955 |
| 7 | 0.3906 | 0.9741 |
| 8 | 0.6504 | 0.9944 |

---

### 🔷 RCWFN — Coulomb Wave Functions
**File:** `src/rcwfn.f`  
**Role:** Computes regular (F_L) and irregular (G_L) Coulomb wave functions for partial wave matching in WAVELJ.

**Physics:**
```
F_L(η, ρ) — regular Coulomb function  (→ sin(ρ - ηln(2ρ) - Lπ/2 + σ_L) as ρ→∞)
G_L(η, ρ) — irregular Coulomb function (→ cos(ρ - ηln(2ρ) - Lπ/2 + σ_L) as ρ→∞)
```

**Method:** Steed's method + continued fractions

**Called by:** WAVELJ (matching), CLINTS (Coulomb integrals)

---

### 🔷 WFGET — Wavefunction Retrieval
**File:** `source.mor` line 32533  
**Role:** Retrieves the precomputed distorted wave u_L(r) at specific grid points by interpolating the stored WAVELJ output.

**Returns:** u_L(r) values at Gauss quadrature points  
**Note:** Returns u_L(r) = r·χ_L(r), not χ_L(r)

---

### 🔷 GETSCT — Compute All S-Matrices
**File:** `source.mor` line 15538  
**Role:** Loops over all partial waves L from 0 to L_max, calls WAVELJ for each L, stores the S-matrix for use in BETCAL.

---

### 🔷 GRDSET — Integration Grid Setup
**File:** `source.mor` line 15710  
**Role:** Sets up the Gauss quadrature grid for the 2D radial DWBA integral over (r_i, r_o). Also sets up the phi_ab angular grid.

**Key variables:**
- `JACOB` — integration Jacobian (r³ for cubic coordinate mapping)
- `RIROWTS = JACOB * ri * ro * w_ij` — combined integration weight

**Coordinate transform:** Uses a cubic mapping r = S³ for efficient sampling, so the Jacobian is JACOB=S³=r.

**Calls:** BSPROD (to evaluate bound state overlaps at grid points), CUBMAP

---

### 🔷 BSPROD — Bound State Product Evaluator
**File:** `source.mor` line 4531  
**Role:** Evaluates products of wavefunctions at integration grid points. Uses Lagrange interpolation (AITLAG) to get wavefunction values between stored grid points.

**ITYPE modes:**
| ITYPE | Computes |
|-------|---------|
| 1 | H(r_i, r_o) = ∫ φ_T(r_t)·V(r_p)·φ_P(r_p) dphi_ab (H-function) |
| 2 | Bound state φ(r) at a point |
| 3 | r·χ(r) at a grid point (scattering wave) |
| 4 | Full scattering wave for tracking |

---

### 🔷 INELDC — Main DWBA Radial Integral
**File:** `source.mor` line 17454  
**Role:** Computes the DWBA transfer T-matrix by numerical double integration over (r_i, r_o).

**Physics (zero-range approximation):**
```
I(Li, Lo, Lx) = ∫∫ u_a*(Li, ri) · H(ri, ro) · u_b(Lo, ro) · RIROWTS · dri · dro

where:
  u_a = r·χ_a(r)  ← distorted wave, incoming channel
  u_b = r·χ_b(r)  ← distorted wave, outgoing channel
  H(ri,ro) = angular kernel integrating over phi_ab
  RIROWTS = JACOB · ri · ro · w_ij  ← integration weight
```

**Full integrand (from source.mor ~line 16548):**
```
TERM = JACOB × ri × ro × weight
Integrand = u_a*(ri) × u_b(ro) × H(ri,ro) × TERM
          = [ri·χ_a(ri)] × [ro·χ_b(ro)] × H × JACOB × ri × ro
          = χ_a(ri) × χ_b(ro) × H × JACOB × ri² × ro²
```

**⚠️ C++ Bug:** C++ stores χ(r) directly (not r·χ(r)), so the `ri²·ro²` term becomes just `ri·ro` — missing one factor of ri and one of ro!

**After integration:** Calls SFROMI to convert I(Li,Lo,Lx) → S-matrix element

---

### 🔷 A12 — Angular Coupling Kernel
**File:** `source.mor` line 1453  
**Role:** Computes the angular momentum transformation coefficient for the DWBA transfer kernel. This is the angular part of the H-function in INELDC.

**Physics:**
```
A12(phi_ab) = XN × xlam(lT,MT) × xlam(lP,MP) 
              × 3j(lT, lP, Lx; 0, 0, 0)
              × xlam(Li, MU) × xlam(Lo, Mx-MU) × √(2Lo+1)
              × 3j(Li, MU, Lo, Mx-MU, Lx, -Mx)
              × cos(MU × phi_ab)

where xlam(l,m) = √((2l+1)/4π) × spherical harmonic factor
```

**Called by:** INELDC (inside the phi_ab angular integration loop)

---

### 🔷 SFROMI — Transfer S-Matrix Assembly
**File:** `source.mor` line 29003  
**Role:** Converts the raw radial integral I(Li,Lo,Lx) into a proper S-matrix element, applying all kinematic and spin factors.

**Physics:**
```
S_sfromi(Li, Lo) = FACTOR × i^(Li+Lo+2Lx+1) × (1/√(2Li+1)) × I(Li,Lo,Lx)
                 × W9J(JBT, 2Lx, JBP, JPI, 2Li, JA, JPO, 2Lo, JB)
                 × √((JPI+1)(JPO+1)(2Lx+1)(JBP+1))

where:
  FACTOR = 2√(ka·kb / (Ecm_a · Ecm_b)) = 0.1110
  JBT = 2×jT = 3 (j=3/2 of neutron in ³⁴Si)
  JBP = 2×jP = 1 (j=1/2 of neutron in deuteron)
  JA = 2 (2×(3/2) for ³³Si ground state J=3/2)
  JB = 1 (2×(1/2) for proton)
```

**Reference value:** |S_sfromi(Li=2, Lo=2, Lx=2)| = 0.01903

**⚠️ C++ Bug:** The 9-J symbol and its statistical factor `√((JPI+1)(JPO+1)(2Lx+1)(JBP+1))` are **missing** from C++ — it only computes the spinless approximation.

---

### 🔷 BETCAL — Beta Amplitude Calculator
**File:** `source.mor` line 3358  
**Role:** Computes the angle-independent β(Lo) amplitudes from S_sfromi elements. These are the partial-wave amplitudes before summing over angles.

**Physics:**
```
β(Lo) = FACTOR_BET × Σ_Li (2Li+1) × CG(Li,0,Lx,Mx,Lo,Mx) 
                    × e^(i(σ_Li + σ_Lo)) × S_sfromi(Li,Lo)

where:
  FACTOR_BET = 0.5/ka = 0.3815
  CG = Clebsch-Gordan coefficient
  σ_L = Coulomb phase = arg(Γ(L+1+iη))
```

**Called by:** XSECTN (before AMPCAL)

---

### 🔷 AMPCAL — Angular Distribution Amplitude
**File:** `source.mor` line 220  
**Role:** Computes the scattering amplitude F(θ) at each angle by summing β(Lo) over partial waves using associated Legendre polynomials.

**Physics:**
```
F(θ, Mx) = Σ_Lo β(Lo, Mx) × P_Lo^Mx(cos θ)
```

**Calls:** PLMSUB (associated Legendre polynomials), EPSLON (Padé extrapolation)

---

### 🔷 XSECTN — Cross Section Calculator
**File:** `source.mor` line 32743  
**Role:** Orchestrates the final cross section computation and output.

**Sequence:**
1. Call BETCAL → get β(Lo) amplitudes
2. Loop over angles → call AMPCAL → get F(θ)
3. Apply XSECTN formula → get dσ/dΩ in mb/sr
4. Call ANAPOW → print output table

**Physics (DWBA cross section formula):**
```
dσ/dΩ = (μa·μb / (ℏ⁴·π)) × (kb/ka) × (1/((2sd+1)(2sA+1))) × Σ_Mx FMNEG × |F(θ,Mx)|²

where:
  μa = reduced mass d+³³Si system
  μb = reduced mass p+³⁴Si system
  sd = deuteron spin = 1     → (2sd+1) = 3
  sA = ³³Si spin = 3/2       → (2sA+1) = 4
  FMNEG = 1 for Mx=0, 2 for Mx>0  (time-reversal degeneracy)
  Conversion: ×10 for fm² → mb
```

---

### 🔷 ANAPOW — Angular Power / Output
**File:** `source.mor` line 578  
**Role:** Formats and prints the final differential cross section dσ/dΩ vs angle table. Also computes analyzing powers if requested.

**Calls:** MUELCO (Mueller matrix coefficients for polarization observables)

---

### 🔷 ELDCS — Elastic Cross Section
**File:** `source.mor` line 12989  
**Role:** Computes the elastic scattering cross section from the optical model S-matrix. Called separately from the transfer calculation.

**Calls:** BETCAL, AMPCAL (same subroutines but for elastic)

---

## 🔗 Data Flow Diagram (Detailed)

```
Input file (test_exact.txt.in)
        │
        ▼
     DATAIN ──────────────────────────────────────────────────┐
     Parses: V,r₀,a,Vi,Vsi,Vso,Rc (both channels)           │
             n,l,j,Ebind (bound states)                       │
             Elab, masses, angles                              │
                                                              │
        │                                                     │
        ├──► BOUND (target: n in ³⁴Si, 0d3/2)               │
        │    │  WOODSX (WS potential on grid)                 │
        │    │  Numerov integration → u(r)                    │
        │    │  depth search for V                            │
        │    └─► stores φ_T(r) = u_T(r)/r                    │
        │                                                     │
        ├──► BOUND (projectile: n in deuteron, 0s1/2)        │
        │    └─► stores φ_P(r) = u_P(r)/r                    │
        │                                                     │
        ├──► WAVSET → GETSCT                                  │
        │         │                                           │
        │         ├──► WAVELJ (L=0..Lmax, incoming)          │
        │         │    RCWFN (Coulomb matching)              │
        │         │    → stores u_a,L(r) = r·χ_a,L(r)       │
        │         │    → S-matrix S_a,L                      │
        │         │                                           │
        │         └──► WAVELJ (L=0..Lmax, outgoing)          │
        │              RCWFN (Coulomb matching)              │
        │              → stores u_b,L(r) = r·χ_b,L(r)       │
        │              → S-matrix S_b,L                      │
        │                                                     │
        └──► GRDSET (set up Gauss quadrature grid)           │
             CUBMAP (cubic coordinate mapping r=S³)           │
             BSPROD (precompute H-functions at grid points)   │
                                                              │
                    ▼                                         │
                 INELDC (main double integral)                │
                 ┌────────────────────────────┐              │
                 │ for each (Li, Lo):          │              │
                 │   WFGET → u_a(Li, ri)      │              │
                 │   WFGET → u_b(Lo, ro)      │              │
                 │   A12   → angular kernel   │              │
                 │   Σ JACOB·ri²·ro²·H·dri·dro│              │
                 │   → raw integral I(Li,Lo,Lx)│              │
                 │   SFROMI → S_sfromi(Li,Lo) │              │
                 └────────────────────────────┘              │
                                                              │
                    ▼                                         │
                 XSECTN                                       │
                    │                                         │
                    ├──► BETCAL                               │
                    │    β(Lo) = FACTOR_BET × Σ_Li            │
                    │           (2Li+1) × CG × e^(iσ)        │
                    │           × S_sfromi(Li,Lo)             │
                    │                                         │
                    ├──► AMPCAL (for each angle θ)           │
                    │    F(θ,Mx) = Σ_Lo β(Lo,Mx)·P_Lo^Mx(cosθ)│
                    │    PLMSUB → Legendre polynomials        │
                    │                                         │
                    └──► ANAPOW                               │
                         dσ/dΩ = prefactor × |F|²            │
                         Output: mb/sr vs degrees             │
```

---

## ⚠️ Known Bugs (C++ Translation)

| # | Bug | Location | Status |
|---|-----|----------|--------|
| 1 | VSI parsing: vsi/rsi0/asi not read | PtolemyParser.cpp | ✅ Fixed |
| 2 | Numerov index offset: wrong u[N] vs u[N-1] | dwba.cpp WavElj() | ✅ Fixed |
| 3 | N_int truncation: was min(NSteps,200) | dwba.cpp | ✅ Fixed |
| 4 | Sign convention: nuclear potential sign | dwba.cpp | ✅ Fixed |
| 5 | SFROMI norm: spurious √((2Li+1)(2Lo+1)) | dwba.cpp | ✅ Fixed |
| 6 | Centrifugal term L(L+1)/r² missing from WavElj | dwba.cpp | ✅ Fixed |
| 7 | **Missing ri²·ro² in INELDC integrand** | dwba.cpp ~line 853 | ❌ Open |
| 8 | **Missing 9-J symbol in SFROMI** | dwba.cpp | ❌ Open |
| 9 | Integral measure: JACOB·ri·ro vs ri²·ro² | dwba.cpp | ❓ Under investigation |
| 10 | Coulomb phase sigma_L in BETCAL | dwba.cpp | ❓ Under investigation |
| 11 | Spectroscopic amplitude SPAM=0.97069 | dwba.cpp | ❓ Under investigation |

---

## 📐 Key Physics Constants (This Reaction)

```
Reaction:    ³³Si(d,p)³⁴Si   Elab = 20 MeV
Q-value:     5.3241 MeV
Ecm_a:       18.843 MeV       (d+³³Si center of mass)
Ecm_b:       24.160 MeV       (p+³⁴Si center of mass)
ka:          1.3109 fm⁻¹      (incoming wavenumber)
kb:          1.0634 fm⁻¹      (outgoing wavenumber)
η_a:         0.704             (Sommerfeld parameter, incoming)
η_b:         0.452             (Sommerfeld parameter, outgoing)
FACTOR:      0.1110            = 2√(ka·kb/(Ecm_a·Ecm_b))
FACTOR_BET:  0.3815            = 0.5/ka
```

---

## 📊 Reference Cross Sections (Ptolemy)

| θ (deg) | dσ/dΩ (mb/sr) |
|---------|--------------|
| 0       | 1.863        |
| 5       | 1.905        |
| 10      | 2.167        |
| 15      | 2.535        |
| 20      | 2.457        |
| 25      | 1.759        |
| 30      | 0.905        |

---

*Last updated: 2026-03-13 — Generated from Ptolemy source.mor analysis*
