# Ptolemy Fortran Source

Original Fortran source code for the Ptolemy DWBA code (April 2007 version), by M.H. Macfarlane and S.C. Pieper, Argonne National Laboratory.

## Quick Start

```bash
# Ubuntu 24.04 — install 32-bit support (recommended)
sudo apt install gcc-multilib gfortran-multilib

# Build
make

# Run
./ptolemy < test01.in > test01.out
```

## Build Requirements

- **gfortran** (GNU Fortran compiler, version 10+)
- **gcc** (for two small C support files)
- **gcc-multilib + gfortran-multilib** (recommended, for 32-bit build)

On Ubuntu 24.04:
```bash
sudo apt install gfortran gcc gcc-multilib gfortran-multilib
```

## Build Modes

### 32-bit build (recommended)

```bash
make ptolemy32
```

The original Ptolemy code was written for 32-bit systems. The internal memory allocator (`GRAB`/`LOC`) stores memory addresses in Fortran `INTEGER` variables (4 bytes). On 64-bit systems, this causes pointer truncation and silent failures. The 32-bit build avoids this entirely.

Requires `gcc-multilib` and `gfortran-multilib` packages.

### 64-bit build

```bash
make ptolemy64
```

Uses a patched `srread_64.c` with `intptr_t` instead of `int` for file descriptor storage. However, the Fortran-side allocator (`gfortran_stuff.f`) still uses `INTEGER` for `LOC()` addresses, which may cause issues with large memory layouts. **Use the 32-bit build if possible.**

### Auto-detect

```bash
make
```

Tries 32-bit first; falls back to 64-bit if multilib packages are not installed.

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
| `gfortran_stuff.f` | Platform-specific routines — GRAB (memory allocator), SECOND (timer), CMPUTR (hostname), LOC |
| `masstable.f` | Nuclear mass table (AME data) |
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
| `srread.c` | Sequential file I/O (original, 32-bit `int` for FILE*) |
| `srread_64.c` | Sequential file I/O (64-bit fix, uses `intptr_t`) |

### Build Infrastructure

| File | Description |
|------|-------------|
| `expand.f` | Macro expander — converts `.mor` files to `.f` (not needed for pre-expanded sources) |
| `macros1` | Macro definitions (selects `com` for fixed common block allocator) |
| `Makefile` | Build system for Ubuntu 24.04 |

## Compilation Flags Explained

```
-std=legacy          Allow old Fortran 77 constructs (GOTO, EQUIVALENCE, etc.)
-fallow-argument-mismatch   Allow type mismatches in subroutine calls
-fno-range-check     Allow integer overflow (used for BOZ constants)
-fallow-invalid-boz  Allow non-standard BOZ literal assignments
-fd-lines-as-comments  Treat 'D' in column 1 as comments (debug lines)
-fno-automatic       Use static storage for local variables (Fortran 77 behavior)
-w                   Suppress all warnings (the code generates many)
```

## Usage

```bash
./ptolemy < input.in > output.out
```

Standard input/output. See `../docs/PTOLEMY_MANUAL.md` for input format and `../docs/PTOLEMY_OUTPUT.md` for output format.

## Known Issues

### 64-bit Pointer Truncation

The `GRAB` subroutine in `gfortran_stuff.f` uses the Fortran `LOC()` intrinsic to get memory addresses, then stores them in `INTEGER` (4-byte) variables. On 64-bit systems, this truncates the address and causes incorrect array indexing. Symptoms: missing output sections, wrong nuclear masses, or silent crashes.

**Fix:** Use the 32-bit build (`make ptolemy32`).

### VSO Factor for Deuterons

Ptolemy divides the spin-orbit coupling by 2S. For S=1 (deuteron), VSO is effectively halved. Published OM parameters for deuterons have this baked in. See `../docs/PTOLEMY_MANUAL.md` §11.1.

### Fortran Carriage Control

Output uses column-1 carriage control characters (`0` = double space, `1` = new page). These appear as literal characters when redirecting to a file. See `../docs/PTOLEMY_OUTPUT.md` §13.

## License

Academic use. Original code by M.H. Macfarlane and S.C. Pieper, Argonne National Laboratory.
