# Ptolemy Fortran Source

Original Fortran source code for the Ptolemy DWBA code (April 2007 version), by M.H. Macfarlane and S.C. Pieper, Argonne National Laboratory.

## Quick Start

```bash
# Ubuntu 24.04
sudo apt install gfortran gcc

# Build and test
make
make test
```

## Build Requirements

- **gfortran** (GNU Fortran compiler, version 10+)
- **gcc** (for two small C helper files)

On Ubuntu 24.04:
```bash
sudo apt install gfortran gcc
```

## Build

```bash
make          # 64-bit build with -O0
make test     # build + run 16O(d,p)17O test
make clean    # remove build artifacts
```

### ⚠️ CRITICAL: Optimization Level

**The code must be compiled with `-O0` on 64-bit systems.** Even `-O1` breaks the output.

The original Ptolemy code relies on Fortran 77 behaviors (implicit SAVE semantics, uninitialized local variables retaining values between calls, pointer-as-integer arithmetic) that modern compilers optimize away at `-O1` and above.

### 32-bit Build (optional)

If you have 32-bit multilib support, the 32-bit build can use `-O2` safely:

```bash
sudo apt install gcc-multilib gfortran-multilib
make ptolemy32
```

The original Ptolemy was developed and tested on 32-bit systems. The pre-built binary in `../digios/analysis/Cleopatra/ptolemy` is 32-bit and statically linked.

## Usage

```bash
./ptolemy < input.in > output.out
```

Standard input/output. See `../docs/PTOLEMY_MANUAL.md` for input format and `../docs/PTOLEMY_OUTPUT.md` for output format.

## Verified Output

The 64-bit `-O0` build produces **identical** results to the reference 32-bit binary:

```
16O(d,p)17O at 20 MeV:  0° DCS = 41.629 mb/sr, TOTAL = 35.050 mb
```

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
