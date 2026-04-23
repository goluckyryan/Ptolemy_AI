# Ptolemy++ — Nuclear Reaction Code (C++)

A C++ translation of the Fortran [Ptolemy](https://www.phy.anl.gov/theory/research/ptolemy/) code for computing distorted-wave Born approximation (DWBA) cross sections for nuclear transfer reactions and optical model elastic scattering.

## Features

### Elastic Scattering ✅
- Spin-0, spin-1/2, and spin-1 projectiles (α, p/n, d/t/³He)
- Woods-Saxon + surface + spin-orbit + Coulomb optical potentials
- Unpolarized differential cross section dσ/dΩ(θ)
- Ratio to Rutherford (σ/σ_Ruth)
- Total reaction cross section
- S-matrix output (magnitude and phase for each L, J)
- Wynn epsilon acceleration (faithful port of Ptolemy's EPSLON subroutine)
- Validated against Fortran Ptolemy (32-bit Cleopatra): **0.08% mean** for both (p,p) and (d,d)

### DWBA Transfer Reactions ✅
- Single-step (d,p), (d,n), (p,d) transfer reactions
- AV18 deuteron bound state wavefunction
- Target bound state with spin-orbit (automatic WS depth search)
- Distorted waves via Numerov integration + Coulomb matching
- GRDSET/InelDc radial integrals (3D grid)
- SFROMI transfer S-matrix with 9J coupling
- XSECTN differential cross sections (CM frame)
- Validated against Fortran Ptolemy: **0.26% mean** (206Hg benchmark, 0–180°)

### Known Limitations
- Single-step transfers only (no coupled channels, no multi-step)
- CM frame output only (no lab-frame Jacobian conversion yet)
- No tensor analyzing powers yet (elastic mode)

## Build

Requires a C++17 compiler. No external dependencies.

### Prerequisites

**Ubuntu/Debian:**
```bash
sudo apt install build-essential gfortran
# For 32-bit Fortran reference binary (optional):
sudo apt install gfortran-multilib gcc-multilib libc6-dev-i386
```

**Fedora/RHEL:**
```bash
sudo dnf install gcc-c++ gcc-gfortran
# For 32-bit: sudo dnf install glibc-devel.i686 libgfortran.i686
```

### C++ Code

```bash
g++ -O2 -std=c++17 -Iinclude \
  src/main.cpp \
  src/dwba/a12.cpp \
  src/dwba/av18_potential.cpp \
  src/dwba/bound.cpp \
  src/dwba/coulin.cpp \
  src/dwba/dwba.cpp \
  src/dwba/grdset_ineldc_faithful.cpp \
  src/dwba/ineldc_collective.cpp \
  src/dwba/math_utils.cpp \
  src/dwba/potential_eval.cpp \
  src/dwba/rcwfn.cpp \
  src/dwba/setup.cpp \
  src/dwba/spline.cpp \
  src/dwba/stubs.cpp \
  src/dwba/thiele_cf.cpp \
  src/dwba/wavelj.cpp \
  src/dwba/xsectn.cpp \
  src/elastic/elastic.cpp \
  src/input/InputGenerator.cpp \
  src/input/Isotope.cpp \
  src/input/Potentials.cpp \
  src/input/PtolemyParser.cpp \
  -o ptolemy_cpp -lm
```

### Fortran Reference Binary (32-bit)

The original Fortran Ptolemy must be built as a **32-bit** binary due to integer/pointer size assumptions in the allocator system.

```bash
cd fortran
make        # builds ./ptolemy (32-bit ELF)
```

Requires `gfortran` with 32-bit multilib support. The `logfac_` size warning during linking is harmless and can be ignored.

⚠️ **Do not build 64-bit** — the code uses `EQUIVALENCE` between `REAL*4` and `REAL*8` arrays with integer packing that breaks on 64-bit.

## Usage

```bash
# From file
./ptolemy++ input.in

# From stdin
./ptolemy++ < input.in

# With angle override
./ptolemy++ input.in 0 180 1
```

Output: DCS table to stdout, diagnostics to stderr.

The parser **auto-detects** elastic vs transfer mode from the reaction line:
- `(d,p)`, `(p,d)`, `(d,n)` → DWBA transfer
- `(d,d)`, `(p,p)`, `(a,a)` → Elastic scattering

## Example Inputs

### Elastic Scattering

```
reset
REACTION: 40Ca(d,d)40Ca(0+ 0.000) ELAB= 20.000
PARAMETERSET dpsb r0target
lmax=30 asymptopia=50

INCOMING
v=90.671 r0=1.150 a=0.762
vi=2.348 ri0=1.334 ai=0.513
vsi=10.218 rsi0=1.378 asi=0.743
vso=3.557 rso0=0.972 aso=1.011
rc0=1.303
;
anglemin=5 anglemax=180 anglestep=5
;
end
```

### DWBA Transfer

```
reset
REACTION: 40Ca(d,p)41Ca(7/2- 0.000) ELAB= 20.000
PARAMETERSET dpsb r0target
lstep=1 lmin=0 lmax=30 maxlextrap=0 asymptopia=50

PROJECTILE
wavefunction av18
r0=1 a=0.5 l=0 rc0=1.2
;
TARGET
JBIGA=0
nodes=0 l=3 jp=7/2
r0=1.25 a=.65
vso=6 rso0=1.10 aso=.65
rc0=1.3
;
INCOMING
v=90.671 r0=1.150 a=0.762
vi=2.348 ri0=1.334 ai=0.513
vsi=10.218 rsi0=1.378 asi=0.743
vso=3.557 rso0=0.972 aso=1.011
rc0=1.303
;
OUTGOING
v=48.457 r0=1.186 a=0.672
vi=2.424 ri0=1.186 ai=0.672
vsi=7.048 rsi0=1.288 asi=0.540
vso=5.284 rso0=0.998 aso=0.590
vsoi=-0.131 rsoi0=0.998 asoi=0.590 rc0=1.349
;
anglemin=0 anglemax=180 anglestep=1
;
end
```

## Input Format

- **Comments:** `$` (inline or line-start) and `'` (line-start)
- **Spacing:** Both `key=value` and `key = value` accepted
- **REACTION line:** `REACTION: Target(proj,eject)Residual(JP Ex) ELAB=energy`
- **Sections:** `PROJECTILE`, `TARGET`, `INCOMING`, `OUTGOING` — each terminated by `;`
- **Binding energy:** Auto-computed from AME2003 masses if not specified

### PARAMETERSET Keywords

| Keyword | Effect |
|---------|--------|
| `dpsb` | High-accuracy grid parameters |
| `r0target` | Use `r0 × A_core^(1/3)` for bound-state radius |
| `tmatch` | 2-point Coulomb matching (Fortran-compatible, default) |
| `wronskian` | 5-point Wronskian matching |
| `lmin=N` / `lmax=N` | Partial wave range |
| `asymptopia=R` | Maximum radius for scattering integration (fm) |

## Project Structure

```
Cpp_AI/
├── src/
│   ├── main.cpp                        # Main driver (elastic / transfer routing)
│   ├── dwba/
│   │   ├── bound.cpp                   # Bound state solver
│   │   ├── wavelj.cpp                  # Distorted wave Numerov integrator
│   │   ├── rcwfn.cpp                   # Coulomb functions (F, G, H±)
│   │   ├── grdset_ineldc_faithful.cpp  # GRDSET + InelDc radial integrals
│   │   ├── xsectn.cpp                  # Transfer cross sections (SFROMI + AMPCAL)
│   │   ├── dwba.cpp                    # DWBA orchestration
│   │   ├── setup.cpp                   # Channel setup and kinematics
│   │   ├── a12.cpp                     # AV18 vertex form factors
│   │   ├── av18_potential.cpp          # AV18 NN potential
│   │   ├── potential_eval.cpp          # Woods-Saxon potential evaluation
│   │   ├── math_utils.cpp              # CG, 6J, 9J, Legendre polynomials
│   │   └── spline.cpp                  # Cubic spline interpolation
│   ├── elastic/
│   │   └── elastic.cpp                 # Elastic scattering solver (Numerov + DCS)
│   └── input/
│       ├── PtolemyParser.cpp           # Ptolemy .in file parser
│       ├── InputGenerator.cpp          # DWBA spec → .in file generator
│       ├── Isotope.cpp                 # Nuclear mass lookup (AME2020)
│       └── Potentials.cpp              # OM potential library (Koning-Delaroche, etc.)
├── include/
│   ├── dwba.h                          # DWBA class
│   ├── elastic.h                       # ElasticSolver class
│   ├── PtolemyParser.h                 # Parser (elastic + transfer)
│   ├── Isotope.h                       # Nuclear isotope data
│   ├── sixj_racah.h / wig9j.h         # Angular momentum coupling
│   └── ...
├── data/
│   ├── mass20.txt                      # AME2020 mass table
│   ├── av18-phi-v                      # AV18 deuteron wavefunction table
│   └── reid-phi-v                      # Reid deuteron wavefunction table
├── docs/                               # Theory documentation
└── fortran_testing/                    # Fortran reference and test scripts
```

## Physics

### Elastic Scattering
Solves the radial Schrödinger equation with complex optical model potential via Numerov integration. S-matrix extracted by matching to Coulomb functions at two asymptotic points. DCS computed from the general spin-dependent scattering amplitude with Clebsch-Gordan coupling.

### DWBA Transfer
1. **Bound states:** Deuteron (AV18) and target (Woods-Saxon + SO) via Numerov with automatic depth search
2. **Distorted waves:** Optical model Numerov integration for incoming and outgoing channels
3. **Radial integrals:** 3D form factor on (r_sum, r_diff, φ) grid with Gauss-Legendre quadrature
4. **Transfer S-matrix:** Angular momentum recoupling via 6J and 9J symbols (SFROMI)
5. **Cross sections:** Coherent partial wave summation (AMPCAL → XSECTN) for dσ/dΩ(θ)

## Validation

| Benchmark | Metric | Result |
|-----------|--------|--------|
| ⁴⁸Ca(p,p) elastic DCS | Mean error vs Cleopatra | **0.08%** |
| ⁴⁸Ca(d,d) elastic DCS | Mean error vs Cleopatra (Vso/2S) | **0.08%** |
| ²⁰⁶Hg(d,p) transfer DCS | Mean error vs Cleopatra (all angles) | **0.26%** |
| ²⁰⁶Hg(d,p) transfer DCS | Peak region (30–90°) | **0.01–0.05%** |

### Fortran Ptolemy Spin-Orbit Convention
Ptolemy has two coupled SO conventions (SDOTL /2S + WOODSX factor 2) that cancel for S=½
but give half-strength SO for S=1 (deuterons). Our C++ uses the correct physics formulas.
To match Ptolemy: divide input Vso by 2S when passing to the solver.
See `docs/PTOLEMY_THEORY.md` and `docs/EPSILON_ALGORITHM.md` for details.

## References

- M.H. Macfarlane and S.C. Pieper, "Ptolemy: A Program for Heavy-Ion Direct-Reaction Calculations," Argonne National Laboratory Report ANL-76-11 (1978)
- R.B. Wiringa, V.G.J. Stoks, and R. Schiavilla, Phys. Rev. C **51**, 38 (1995) — AV18 potential

## License

Academic use. Original Ptolemy code by M.H. Macfarlane and S.C. Pieper (Argonne National Laboratory).
