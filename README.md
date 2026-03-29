# Ptolemy C++ — DWBA Transfer Reaction Code

A C++ translation of the Fortran [Ptolemy](https://www.phy.anl.gov/theory/research/ptolemy/) code for computing distorted-wave Born approximation (DWBA) cross sections for nuclear transfer reactions.

## Status

**Work in progress.** Single-step (d,p) transfer reactions with AV18 deuteron wavefunctions are functional. Accuracy vs the Fortran Cleopatra binary is typically **1–5%** at forward angles, with larger deviations at diffraction minima and extreme back angles.

### What works
- Input parser: standard Ptolemy `.in` format (including `$` comments, spaced `key = value`)
- AV18 deuteron bound state (automatic depth search)
- Target bound state with SO potential (automatic WS depth search)
- Optical model scattering (Woods-Saxon + volume/surface imaginary + spin-orbit)
- Distorted waves via Numerov integration + Coulomb matching
- GRDSET/InelDc radial integrals (3D grid: sum, difference, phi coordinates)
- SFROMI: transfer S-matrix with 9J coupling
- XSECTN: differential cross sections (CM frame)
- Two scattering-wave matching methods: Wronskian (default) and TMATCH (Fortran-compatible)

### Known limitations
- Single-step transfers only (no coupled channels, no multi-step)
- CM frame output only (no lab-frame Jacobian conversion yet)
- No elastic cross section output
- ~2–8% DCS error vs Fortran, mainly from radial integral differences at intermediate partial waves
- Large-L (>28) integrals have sub-barrier truncation errors (physically negligible: <0.1% of DCS)

## Build

Requires a C++17 compiler. No external dependencies.

```bash
g++ -O2 -std=c++17 -Iinclude \
  main.cpp \
  src/input/PtolemyParser.cpp \
  src/dwba/grdset_ineldc_faithful.cpp \
  src/dwba/wavelj.cpp \
  src/dwba/rcwfn.cpp \
  src/dwba/math_utils.cpp \
  src/dwba/bound.cpp \
  src/dwba/dwba.cpp \
  src/dwba/setup.cpp \
  src/dwba/xsectn.cpp \
  src/dwba/a12.cpp \
  src/dwba/av18_potential.cpp \
  src/dwba/potential_eval.cpp \
  src/dwba/spline.cpp \
  src/input/Isotope.cpp \
  -o ptolemy_cpp -lm
```

## Usage

```bash
./ptolemy_cpp < input.in
```

Output goes to stdout (DCS table), diagnostics to stderr.

### Example input

```
$============================================ Ex=0.87(1s1/2)AK
reset
REACTION: 16O(d,p)17O(1/2+ 0.87) ELAB= 20.000
PARAMETERSET dpsb r0target
lstep=1 lmin=0 lmax=30 maxlextrap=0 asymptopia=50

PROJECTILE
wavefunction av18
r0=1 a=0.5 l=0 rc0=1.2
;
TARGET
JBIGA=0
nodes=1 l=0 jp=1/2
r0=1.25 a=.65
vso=6 rso0=1.10 aso=.65
rc0=1.3
;
INCOMING
v = 88.955 r0 = 1.149 a = 0.751
vi = 2.348 ri0 = 1.345 ai = 0.603
vsi = 10.218 rsi0 = 1.394 asi = 0.687
vso = 3.557 rso0 = 0.972 aso = 1.011
rc0 = 1.303
;
OUTGOING
v = 49.870 r0 = 1.146 a = 0.675
vi = 1.959 ri0 = 1.146 ai = 0.675
vsi = 7.758 rsi0 = 1.302 asi = 0.528
vso = 5.314 rso0 = 0.934 aso = 0.590
vsoi = -0.100 rsoi0 = 0.934 asoi = 0.590
rc0 = 1.419
;
anglemin=0 anglemax=180 anglestep=1
;
end
```

### Input format notes

- **Comments:** `$` (inline or line-start) and `'` (line-start) are comment characters
- **Spacing:** Both `key=value` and `key = value` are accepted
- **REACTION line:** `REACTION: Target(proj,eject)Residual(JP Ex) ELAB=energy`
- **Sections:** `PROJECTILE`, `TARGET`, `INCOMING`, `OUTGOING` — each terminated by `;`
- **PARAMETERSET:** `dpsb` is the standard high-accuracy parameter set
- **Angle output:** Always CM frame. Set `anglemin`, `anglemax`, `anglestep`.
- **Binding energy:** Auto-computed from nuclear masses if not specified

### Keywords on PARAMETERSET line

| Keyword | Effect |
|---------|--------|
| `dpsb` | High-accuracy grid parameters (NPSUM=40, NPDIF=40, NPPHI=20) |
| `r0target` | Use `r0 * A_core^(1/3)` for bound-state radius |
| `tmatch` | 2-point Coulomb matching (Fortran-compatible) |
| `wronskian` | 5-point Wronskian matching (default) |
| `lstep=N` | Partial wave step |
| `lmin=N` / `lmax=N` | Partial wave range |
| `asymptopia=R` | Maximum radius for scattering integration (fm) |
| `maxlextrap=N` | Max L for extrapolation (0 = no extrapolation) |

### Matching methods

- **Wronskian** (default): 5-point stencil at `NSTEP-3`. Better for some cases (e.g., L=0 transfer at 10 MeV: 0.6% mean error).
- **TMATCH**: 2-point matching at `NSTEP` and `NSTEP-NBACK`. Matches Fortran's SUMMAX exactly. Better for L=2 ground-state transfers (~2.8% mean error).

Choose based on which gives better agreement for your specific reaction.

## Project structure

```
Cpp_AI/
├── main.cpp                  # Main driver (parses input, runs DWBA)
├── include/
│   ├── dwba.h                # DWBA class definition
│   ├── Isotope.h             # Nuclear mass table
│   ├── PtolemyParser.h       # Input file parser
│   ├── av18.h / av18_potential.h  # AV18 deuteron wavefunction
│   ├── sixj_racah.h          # 6J / Racah coefficients
│   ├── wig9j.h               # 9J symbols
│   └── ...
├── src/
│   ├── dwba/
│   │   ├── bound.cpp         # Bound state solver (WS depth search + Numerov)
│   │   ├── wavelj.cpp        # Distorted wave solver (Numerov + Coulomb matching)
│   │   ├── rcwfn.cpp         # Regular/irregular Coulomb functions
│   │   ├── grdset_ineldc_faithful.cpp  # GRDSET + InelDc radial integrals
│   │   ├── xsectn.cpp        # Cross section computation (SFROMI + AMPCAL)
│   │   ├── dwba.cpp          # DWBA orchestration
│   │   ├── setup.cpp         # Channel setup and kinematics
│   │   ├── a12.cpp           # AV18 vertex form factors
│   │   ├── av18_potential.cpp # AV18 potential evaluation
│   │   ├── potential_eval.cpp # Woods-Saxon potential evaluation
│   │   ├── math_utils.cpp    # Clebsch-Gordan, 6J, 9J, Legendre
│   │   └── spline.cpp        # Cubic spline interpolation
│   └── input/
│       ├── PtolemyParser.cpp  # Ptolemy .in file parser
│       └── Isotope.cpp        # Nuclear mass lookup (AME2020)
├── data/
│   └── mass20.txt            # AME2020 mass table
├── docs/                     # Theory documentation
├── fortran_testing/           # Fortran reference code and test scripts
│   ├── source_annotated.f    # Annotated Ptolemy Fortran source
│   ├── ptolemy_annotated      # Annotated Fortran binary (with debug prints)
│   └── build_annotated.sh    # Build script for annotated binary
└── tests/                    # Unit and comparison tests
```

## Physics overview

The code computes DWBA differential cross sections for direct nuclear transfer reactions. For a (d,p) reaction:

1. **Bound states:** Solve for the deuteron (AV18 NN potential) and target (Woods-Saxon) bound-state wavefunctions via Numerov integration with automatic well-depth search.

2. **Distorted waves:** Solve the Schrödinger equation with optical model potentials for incoming (deuteron) and outgoing (proton) channels. Coulomb + nuclear potential with spin-orbit coupling.

3. **Radial integrals:** Compute the 3D form factor integral over (r_in, r_out, φ) using the GRDSET grid with Gauss-Legendre quadrature and spline interpolation.

4. **Transfer S-matrix:** Apply angular momentum coupling (6J, 9J symbols) via SFROMI to convert radial integrals into S-matrix elements labeled by (JP, JT, LX, L_out−L_in).

5. **Cross sections:** Sum coherently over partial waves with Legendre polynomials (AMPCAL) to produce dσ/dΩ(θ) in the CM frame.

## Fortran reference

The Fortran Cleopatra binary is at:
```
digios/analysis/Cleopatra/ptolemy
```
This is the 32-bit clean Ptolemy build used as the reference for validation.

## License

Academic use. Original Ptolemy code by M.H. Macfarlane and S.C. Pieper (Argonne National Laboratory).
