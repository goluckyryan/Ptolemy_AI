# Ptolemy Fortran Source

Original Fortran source code for the Ptolemy DWBA code (April 2007 version), by M.H. Macfarlane and S.C. Pieper, Argonne National Laboratory.

## ⚠️ WARNING: 64-bit Builds Are Not Reliable

**Although this code compiles and runs on 64-bit systems, the 64-bit build can produce incorrect results.** The original code was developed for 32-bit platforms and relies on Fortran 77 behaviors (implicit SAVE, integer pointer arithmetic, uninitialized variable persistence) that do not translate reliably to 64-bit, even at `-O0`.

**Always use the 32-bit build for production work.** The 64-bit build (`make ptolemy64`) is provided only for adding debug prints and tracing internal variables.

For a pre-built, verified 32-bit binary, use `../digios/analysis/Cleopatra/ptolemy`.

## Quick Start

```bash
# Ubuntu 24.04 — install 32-bit multilib support
sudo apt install gfortran gcc gcc-multilib gfortran-multilib

# Build (32-bit, recommended) and test
make
make test
```

## Build Requirements

- **gfortran** (GNU Fortran compiler, version 10+)
- **gcc** (for two small C helper files)
- **gcc-multilib / gfortran-multilib** (for 32-bit build, recommended)

On Ubuntu 24.04:
```bash
sudo apt install gfortran gcc gcc-multilib gfortran-multilib
```

## Build

```bash
make             # 32-bit build (RECOMMENDED)
make ptolemy64   # 64-bit build (DEBUG ONLY — results may be incorrect)
make test        # build + run 16O(d,p)17O test
make clean       # remove build artifacts
```

### 32-bit Build (default, recommended)

The default `make` target builds a 32-bit binary with `-O2`, matching the original development platform. Requires `gcc-multilib` and `gfortran-multilib`.

### 64-bit Build (debugging only)

`make ptolemy64` builds a 64-bit binary with `-O0`. This is useful for adding `PRINT` statements to trace internal variables, but **should never be used as a reference for cross section values**. Even at `-O0`, some reactions produce different results from the 32-bit build.

## Usage

```bash
./ptolemy < input.in > output.out
```

Standard input/output. See `../docs/PTOLEMY_MANUAL.md` for input format and `../docs/PTOLEMY_OUTPUT.md` for output format.

## Verified Output

The **32-bit** build produces identical results to the reference pre-built binary:

```
16O(d,p)17O at 20 MeV:  0° DCS = 41.629 mb/sr, TOTAL = 35.050 mb
```

The 64-bit build may reproduce this for simple test cases, but can diverge for other reactions. Do not rely on it as a reference.

## Source Files

### Core Ptolemy

| File | Description |
|------|-------------|
| `ptolemy-main.f` | Main program — calls CONTRL, orchestrates the calculation |
| `source.f` | **Main source** (~38k lines) — CONTRL, DATAIN, BOUND, WAVELJ, GRDSET, INELDC, SFROMI, XSECTN, and all subroutines |
| `fortlib.f` | Math library — Gauss-Legendre quadrature, Legendre polynomials, Racah/6J/9J symbols, splines |

### Support Libraries

| File | Description |
|------|-------------|
| `gfortran_stuff.f` | Platform-specific routines — GRAB (memory allocator), SECOND (timer), CMPUTR (hostname) |
| `masstable.f` | Nuclear mass table (AME data, compiled into DATA statements) |
| `numbered_store.f` | Named array storage manager (NALLOC/NFREE) |
| `linkule.f` | Bound state wavefunction plugins (AV18, Reid, Woods-Saxon, etc.) |
| `linkulesfitters.f` | Optical model parameter libraries (global OM fits) |
| `av18.f` | Argonne v18 NN potential — deuteron wavefunction |
| `keep.f` / `keepsub.f` / `keptsub.f` | KEEP/NSDUMP — save/restore numbered storage |
| `phiffer.f` | PHIFFER — phi-dependent form factor code |

### C Support

| File | Description |
|------|-------------|
| `dtime.c` | CPU timer (calls `getrusage`) |
| `srread.c` | Sequential file I/O for KEEP/NSDUMP |

### Build Infrastructure

| File | Description |
|------|-------------|
| `expand.f` | Macro expander — converts `.mor` source to `.f` (not needed; pre-expanded `.f` files are provided) |
| `macros1` | Macro definitions used by `expand.f` (selects `com` for fixed common block allocator) |

## Compilation Flags

```
-O0                      REQUIRED on 64-bit (code breaks at -O1+)
-std=legacy              Allow old F77 constructs (assigned GOTO, Hollerith, etc.)
-fallow-argument-mismatch  Allow type mismatches in subroutine calls
-fno-range-check         Allow integer overflow (BOZ constants)
-fallow-invalid-boz      Allow non-standard BOZ literal assignments
-fd-lines-as-comments    Treat 'D' in column 1 as comments
-w                       Suppress warnings (the code generates hundreds)
```

## Known Issues

### Optimization Breaks Output

At `-O1` and above, the compiler optimizes away implicit SAVE behavior and/or reorders operations involving integer/pointer arithmetic. This causes the mass table lookup to fail silently (returning Z=0 for some nuclei) and the DCS calculation to never execute.

### VSO Factor for Deuterons

Ptolemy divides the spin-orbit coupling by 2S. For S=1 (deuteron), VSO is effectively halved. Published OM parameters (e.g., An & Cai 2006) have this baked in. See `../docs/PTOLEMY_MANUAL.md` §11.1.

### Fortran Carriage Control

Output uses column-1 carriage control characters (`0` = double space, `1` = new page, `+` = overprint). These appear as literal characters when redirecting to a file. See `../docs/PTOLEMY_OUTPUT.md` §13.

### 64-bit Pointer Arithmetic

The `GRAB` subroutine uses `LOC()` to compute array offsets. The `gfortran_stuff.f` in this directory has `INTEGER*8` declarations for the `LOC()` variables, ensuring correct pointer arithmetic on 64-bit systems. The original code used default `INTEGER` (4 bytes).

## License

Academic use. Original code by M.H. Macfarlane and S.C. Pieper, Argonne National Laboratory.
