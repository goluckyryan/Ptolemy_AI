# Ptolemy User Manual

> **Version:** Ptolemy (Argonne National Laboratory, Fortran)  
> **Purpose:** DWBA and coupled-channels calculation of direct nuclear reaction cross sections  
> **Authors:** R. Schiff, S. C. Pieper, R. P. Goddard (ANL)  
> **This document:** Practical user guide compiled from source code, input examples, and validation work at ANL (2025–2026)

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Running Ptolemy](#2-running-ptolemy)
3. [Input File Structure](#3-input-file-structure)
4. [Global Control Keywords](#4-global-control-keywords)
5. [Reaction Definition — Old-Style (CHANNEL)](#5-reaction-definition--old-style-channel)
6. [Reaction Definition — New-Style (REACTION)](#6-reaction-definition--new-style-reaction)
7. [Potential Parameters](#7-potential-parameters)
8. [Radius Conventions](#8-radius-conventions)
9. [Bound State Options](#9-bound-state-options)
10. [Output Interpretation](#10-output-interpretation)
11. [Known Quirks and Bugs](#11-known-quirks-and-bugs)
12. [Complete Input Examples](#12-complete-input-examples)
13. [Quick Reference Card](#13-quick-reference-card)

---

## 1. Introduction

Ptolemy computes differential cross sections for direct nuclear reactions using the **Distorted Wave Born Approximation (DWBA)**. It supports:

- **Elastic scattering**: (p,p), (d,d), (α,α), (³He,³He), etc.
- **Transfer reactions**: (d,p), (p,d), (³He,d), (t,p), etc.
- **Coupled-channels** (inelastic): limited support

The calculation proceeds as:
1. Solve the optical model Schrödinger equation in the entrance and exit channels → distorted waves χ_a(r), χ_b(r)
2. Compute bound state wavefunctions φ_T(r) (target) and φ_P(r) (projectile)
3. Evaluate the DWBA transfer amplitude (finite-range or zero-range)
4. Sum partial waves → dσ/dΩ as function of CM angle

**Key physics inputs:**
- Optical Model Potentials (OMP): Woods-Saxon central + surface + spin-orbit + Coulomb
- Bound state: Woods-Saxon potential adjusted to give correct binding energy
- Spectroscopic factors (from shell model or experiment)

---

## 2. Running Ptolemy

```bash
./ptolemy < input.in > output.out
```

Ptolemy reads from **stdin** and writes to **stdout**. Redirect both.

**Typical workflow:**
```bash
# Run calculation
./ptolemy < o16dp_gs.in > o16dp_gs.out

# Check for errors
grep -i "error\|warning\|****" o16dp_gs.out

# Extract cross section
grep -A 200 "CROSS SECTION" o16dp_gs.out | head -50
```

---

## 3. Input File Structure

### 3.1 General Format

Ptolemy input is **free-format keyword-driven**. Each keyword is on its own line (or `KEYWORD = VALUE`). Keywords are **case-insensitive** and can be abbreviated to 8 characters.

```
KEYWORD
KEYWORD = VALUE
KEYWORD VALUE
```

### 3.2 Comments

Lines beginning with `$` are comments:
```
$ This is a comment
HEADER  60Ni(p,p) at 30 MeV   $ HEADER is not a comment line
```

Lines beginning with `0INPUT...` in the output echo the input.

### 3.3 Section Delimiters

| Symbol | Meaning |
|--------|---------|
| `;`    | End of a block (channel definition, bound state, etc.) |
| `END`  | End of entire input — Ptolemy stops |
| `RETURN` | Same as END |
| `RESET` | Reset all parameters to defaults (start fresh) |
| `QUIT` | Abort immediately |

### 3.4 Two Input Styles

**Old-style** (classic Ptolemy):
```
CHANNEL  2H + 60Ni
ELAB = 60.0
V = 81.919  R0 = 1.150  A = 0.768
...
ELASTIC SCATTERING
;
END
```

**New-style** (used by Cleopatra, more structured):
```
reset
REACTION: 16O(d,p)17O(5/2+ 0.000) ELAB=20.0
PARAMETERSET dpsb ...
INCOMING
  V=88.955  R0=1.149  A=0.751
  ...
;
OUTGOING
  ...
;
TARGET
  nodes=0 l=2 jp=5/2
  v=60.0 r0=1.10 a=0.65
;
end
```

---

## 4. Global Control Keywords

### 4.1 Header and Labels

| Keyword | Alias | Description |
|---------|-------|-------------|
| `HEADER text` | `UEBERSCH` | Label printed in output header. Up to ~70 chars. |
| `TITLE text` | — | Alternative label |
| `NEWPAGE` | — | Force page break in output |

### 4.2 Flow Control

| Keyword | Description |
|---------|-------------|
| `END` | Stop processing; exit |
| `RETURN` | Same as END |
| `RESET` | Reset all parameters to defaults |
| `QUIT` | Abort immediately |
| `LISTKEYS` | Print all available keywords to output |

### 4.3 Angular Range

| Keyword | Alias | Description | Default |
|---------|-------|-------------|---------|
| `ANGLEMIN = value` | `WINKELMI` | Minimum lab angle (degrees) | 0 |
| `ANGLEMAX = value` | `WINKELMA` | Maximum lab angle (degrees) | 180 |
| `ANGLESTEP = value` | `SCHRITTE` | Step in lab angle (degrees) | 1 |
| `LABANGLES` | — | Output in lab frame (new-style) | CM |
| `CMANGLES` | — | Output in CM frame | default |

In new-style input, can also specify inline:
```
anglemin=0 anglemax=180 anglestep=5
```

### 4.4 Numerical Accuracy

| Keyword | Alias | Description | Default |
|---------|-------|-------------|---------|
| `ACCURACY = value` | `GENAUIGK` | Convergence tolerance for bound state search | 1e-5 |
| `STEPSPER = value` | `SCHRITTL` | Steps per asymptotic range (bound state grid) | 8 |
| `ASYMPTOPIA = value` | `ASYMPTOP` | Maximum r for integration (fm) | auto |
| `STEPSIZE = value` | `SCHRITTL` | Fixed step size (fm); overrides STEPSPER | auto |
| `WRITESTEP = value` | `SCHREIBS` | Print bound state wf every N fm (0 = no print) | 0 |
| `LMAX = value` | `LMAXHINZ` | Maximum L for partial wave sum | auto |
| `LMIN = value` | `LMINABZI` | Minimum L | 0 |
| `MAXLEXTRAP = value` | — | Extrapolation beyond LMAX | 0 |

### 4.5 Output Options

| Keyword | Description |
|---------|-------------|
| `PRINT = 1001` | Print elastic S-matrix in JP basis (see §10.3) |
| `CROSSSEC` | Compute and print cross section (new-style) |
| `WRITENS` | Write NS (normalization) table |
| `LABANGLE` | Output lab-frame angles |
| `CMANGLES` | Output CM-frame angles (default) |
| `PRINTWAVE` | After BOUNDSTATE: print potential + wavefunction table |
| `COMPLEXW` | Use complex W (volume imaginary) instead of surface W |
| `BATCH` | Suppress some interactive prompts |

---

## 5. Reaction Definition — Old-Style (CHANNEL)

### 5.1 CHANNEL Keyword

```
CHANNEL  <A>H + <Z>X
```

Defines projectile + target. Uses standard nuclear notation:
- `1H` = proton, `2H` = deuteron, `3H` = triton
- `3He`, `4He` = helium isotopes
- `16O`, `60Ni`, `148Sm`, etc.

Examples:
```
CHANNEL  1H + 60Ni      $ proton on 60Ni
CHANNEL  2H + 16O       $ deuteron on 16O
CHANNEL  4He + 148Sm    $ alpha on 148Sm
```

### 5.2 Energy

```
ELAB = value    $ Lab-frame kinetic energy (MeV)
```

### 5.3 Elastic Scattering Block

```
CHANNEL  1H + 60Ni
ELAB = 30.0
[potential parameters]
ELASTIC SCATTERING    $ or just ELASTIC
;
```

The semicolon `;` ends the channel block and triggers computation.

### 5.4 Transfer Reaction Block (Old-Style)

```
CHANNEL  2H + 16O       $ entrance channel
ELAB = 20.0
r0target
[incoming OMP parameters]
TRANSFER
JBIGA = 0+              $ spin-parity of target A (16O: 0+)

CHANNEL  1H + 17O       $ exit channel
ELAB = 0.0              $ Ptolemy computes from Q value
[outgoing OMP parameters]
JBIGB = 5/2+            $ spin-parity of residual B (17O g.s.)

BOUNDSTATE  0  2  5  4.143    $ nodes l 2j BE(MeV)
[bound state potential]
;
END
```

### 5.5 JBIGA / JBIGB

```
JBIGA = 0+      $ J^pi of target nucleus (entrance channel)
JBIGB = 5/2+    $ J^pi of residual nucleus (exit channel)
```

Format: `J pi` where pi is `+` or `-`. Examples: `0+`, `1/2+`, `3/2-`, `5/2+`

---

## 6. Reaction Definition — New-Style (REACTION)

Used by **Cleopatra** and preferred for transfer reactions.

### 6.1 REACTION Line

```
REACTION: 16O(d,p)17O(5/2+ 0.000) ELAB=20.0
```

Format: `Target(projectile,ejectile)Residual(Jpi Ex) ELAB=value`
- Ex = excitation energy in MeV
- Jpi = spin-parity of residual state

### 6.2 PARAMETERSET

Controls global calculation options:
```
PARAMETERSET dpsb labangles r0target lstep=1 lmin=0 lmax=40 maxlextrap=0 asymptopia=30
```

Common options:
| Option | Description |
|--------|-------------|
| `dpsb` | Standard (d,p) parameter set |
| `r0target` | Use target-only radius convention |
| `labangles` | Output in lab frame |
| `lmin=N` | Minimum partial wave |
| `lmax=N` | Maximum partial wave |
| `asymptopia=N` | Max integration radius (fm) |

### 6.3 INCOMING / OUTGOING Blocks

Optical model potentials for entrance and exit channels:
```
INCOMING
V=88.955   R0=1.149  A=0.751
VI=2.348   RI0=1.345 AI=0.603
VSI=10.218 RSI0=1.394 ASI=0.687
VSO=7.114  RSO0=0.972 ASO=1.011
RC0=1.303
;
```

### 6.4 TARGET Block

Bound state of transferred particle in residual nucleus:
```
TARGET
nodes=0 l=2 jp=5/2
v=60.0 r0=1.10 a=0.65
rc0=1.30
writestep=0.1
;
```

Parameters: `nodes` (quantum number n), `l`, `jp` (= 2j), binding energy determined from REACTION line.

### 6.5 PROJECTILE Block

Bound state of transferred particle in projectile (for finite-range):
```
PROJECTILE
wavefunction av18 r0=1 a=0.5 l=0
;
```

Options:
- `wavefunction av18` — use AV18 potential for deuteron internal wf
- `wavefunction file` — read from external file
- `r0=N a=N l=N` — simple WS parameters

---

## 7. Potential Parameters

All potentials are **Woods-Saxon type**. The physical radius is computed as `R = r0 × A^(1/3)` (with convention set by CHANNEL or PARAMETERSET).

### 7.1 Complete Parameter Table

| Keyword | Alias | Description | Units |
|---------|-------|-------------|-------|
| `V` | — | Real central WS depth | MeV |
| `R0` | — | Real central radius parameter | fm |
| `A` | — | Real central diffuseness | fm |
| `VI` | — | Imaginary volume depth | MeV |
| `RI0` | — | Imaginary volume radius | fm |
| `AI` | — | Imaginary volume diffuseness | fm |
| `W` | — | Imaginary volume depth (alias for VI) | MeV |
| `VSI` | — | Imaginary surface depth (derivative WS) | MeV |
| `RSI0` | `RSI` | Imaginary surface radius | fm |
| `ASI` | — | Imaginary surface diffuseness | fm |
| `VSO` | — | Real spin-orbit depth | MeV |
| `RSO0` | `RSO` | Real spin-orbit radius | fm |
| `ASO` | — | Real spin-orbit diffuseness | fm |
| `VSOI` | — | Imaginary spin-orbit depth | MeV |
| `RSOI0` | `RSOI` | Imaginary spin-orbit radius | fm |
| `ASOI` | — | Imaginary spin-orbit diffuseness | fm |
| `RC0` | `RC` | Coulomb radius parameter | fm |

### 7.2 Potential Forms

**Volume central** (type 1):
```
V_WS(r) = -V / (1 + exp((r - R)/a))
```

**Imaginary volume** (`COMPLEXW` or `VI`):
```
W_vol(r) = -VI / (1 + exp((r - RI)/aI))
```

**Imaginary surface** (derivative WS, `VSI`):
```
W_surf(r) = -4 VSI × a_SI × d/dr [1/(1 + exp((r - RSI)/aSI))]
           = VSI × exp((r-RSI)/aSI) / (1 + exp((r-RSI)/aSI))²
```
(Surface form — peaks at r = RSI)

**Spin-orbit** (`VSO`, type 2):
```
V_SO(r) = VSO × (2/(ASO × r)) × d/dr [1/(1 + exp((r - RSO)/aSO))] × (L·S)
         = -VSO × (2/ASO) × exp/((1+exp)² × r)  × (L·S)
```

**Coulomb** (uniform sphere):
```
V_C(r) = Z₁Z₂e²/(2RC) × (3 - r²/RC²)   r < RC
V_C(r) = Z₁Z₂e²/r                        r > RC
```
where `RC = rc0 × A_target^(1/3)` (r0target convention).

### 7.3 Notes on VI vs W

- `W` and `VI` are synonyms for imaginary volume depth
- If both `VI` and `VSI` are nonzero and have **different** R0/a values, use **two separate potential calls** in C++:
  ```cpp
  s.AddVolumeWS({V, 0}, R0_V, a_V);    // real part only
  s.AddVolumeWS({0, -VI}, R0_I, a_I);  // imaginary part only
  ```

---

## 8. Radius Conventions

Ptolemy supports several conventions for converting the radius parameter `r0` to a physical radius `R`:

| Keyword | Formula | When used |
|---------|---------|-----------|
| `r0target` | `R = r0 × A_target^{1/3}` | Cleopatra default; light projectiles |
| `R0DEFAULT` (default) | `R = r0 × A_t^{1/3}` for A_proj ≤ 2.5; `r0 × (A_t^{1/3} + A_p^{1/3})` for A_proj > 2.5 | Ptolemy internal default |
| `R0SUM` | `R = r0 × (A_t^{1/3} + A_p^{1/3})` | Both masses always |
| `R0MTOT` | `R = r0 × (A_t + A_p)^{1/3}` | Total mass |

**Practical rule:**
- For (p,p), (d,p), (d,d): `r0target` and `R0DEFAULT` give the **same result** (A_proj ≤ 2)
- For (α,α), (³He,d): they **differ** — `r0target` gives smaller radius
- **Cleopatra always writes `r0target`** → always use `r0target` when matching Cleopatra output

---

## 9. Bound State Options

### 9.1 Old-Style BOUNDSTATE

```
BOUNDSTATE  n  l  jp  BE
[potential block]
;
```

Parameters:
| Parameter | Description | Example |
|-----------|-------------|---------|
| `n` | Number of radial nodes (0,1,2,...) | `0` |
| `l` | Orbital angular momentum | `2` (d-wave) |
| `jp` | Twice the total angular momentum (2j) | `5` (j=5/2) |
| `BE` | Binding energy (MeV, positive) | `4.143` |

The potential depth `V` is adjusted by bisection to reproduce the binding energy exactly.

Example — 17O g.s. (0d₅/₂, BE=4.143 MeV):
```
BOUNDSTATE  0  2  5  4.143
V=60.0   R0=1.10  A=0.65
RC0=1.30
;
```

### 9.2 New-Style TARGET/PROJECTILE

```
TARGET
nodes=0 l=2 jp=5/2
v=60.0 r0=1.10 a=0.65
rc0=1.30
;
```

Note: `jp=5/2` means j=5/2 (not doubled); `jp=2` means j=2 (integer j); internally Ptolemy always works with doubled integers.

### 9.3 Wavefunction Output

Add `writestep = 0.1` to print φ(r) = u(r)/r every 0.1 fm in the output:

```
BOUNDSTATE  0  2  5  4.143
V=60.0  R0=1.10  A=0.65
writestep = 0.1
;
```

Output format:
```
 PSI=F(R)/R
    STEP      RADIUS        WAVEFUNCTION
       1     0.00000     0.0000
       2     0.08125     0.1077E-02
       ...
```

Use `PRINTWAVE` after a `BOUNDSTATE` block (with `writestep` set) to also print the effective potential.

### 9.4 Grid Parameters

Ptolemy automatically determines the grid:
- Step size: `h = min(1/κ, a) / STEPSPER` where κ = sqrt(2μBE)/ℏ
- Max r: controlled by `ASYMPTOPIA` (default ~30 fm, auto-extended if needed)
- STEPSPER default = 8 (8 steps per asymptotic range)

---

## 10. Output Interpretation

### 10.1 Cross Section Table

The main output table for both elastic and transfer:

```
  LAB.  REACTION     REACTION   LOW L  HIGH L   % FROM
                                                          INCOMING ELASTIC      OUTGOING    REACTION ...
  ANGLE  LAB. MB        /RUTH     %/L   % ERROR  L>LMAX
```

**Column mapping** (0-based):
| Col | Header | Meaning |
|-----|--------|---------|
| 0 | C.M. | CM angle (degrees) |
| 1 | LAB | Lab angle (degrees) |
| 2 | RUTHERFORD | σ/σ_Rutherford |
| 3 | C.M. (σ) | **σ_CM (mb/sr) ← use this** |
| 4 | LAB (σ) | σ_LAB (mb/sr) |
| 5 | C.M. (σ_Ruth) | σ_Rutherford CM |

> ⚠️ **Always extract column 3 (σ_CM) for comparison with CM-frame theory codes.**  
> Column 4 (σ_LAB) is the lab-frame value and differs at backward angles.

### 10.2 Elastic Scattering Output

```
  0.00   41.629      0.000000   14.15    0.00     0.0   ...
  5.00   40.190      0.000842   13.42    0.00     0.0   ...
```

For elastic, `REACTION LAB. MB` column = σ_CM.

### 10.3 S-Matrix Output (PRINT=1001)

Add `PRINT=1001` to get S-matrix in JP (j=l±½) basis:

```
ELASTIC S-MATRIX FOR L = 3, JP = 5/2: Re +- Im I
```

Format: one line per (L, J) pair:
- `JP = 2L-1` → J = L-½
- `JP = 2L+1` → J = L+½

> ⚠️ The **normal** Ptolemy output table (L L' LX columns) is in the **LX spin-flip basis**, NOT the JP basis. Use `PRINT=1001` for the real S-matrix elements.

### 10.4 Bound State Output

```
0        TARGET BOUND STATE PARAMETERS
0E =   -4.1431 MEV     KAPPA = 0.43368
```

Ptolemy prints:
- Final V found (depth in MeV)
- Binding energy E (negative)
- κ = sqrt(2μBE)/ℏ (inverse decay length, fm⁻¹)
- Wavefunction table if `writestep` is set

### 10.5 Common Output Messages

| Message | Meaning |
|---------|---------|
| `**** WARNING: SEARCH FOR DIF GRIDS STOPPED` | Grid extended beyond ASYMPTOPIA; usually OK |
| `**** ERROR: N BOUND STATE WAVEFUNCTIONS NEEDED BEYOND ASYMPTOPIA` | Increase ASYMPTOPIA |
| `NO SOLUTION FOUND` | Bisection failed; check potential parameters |
| `WAVEFUNCTION HAS N NODES BUT M ARE DESIRED` | Wrong node count; adjust V range |

---

## 11. Known Quirks and Bugs

### 11.1 VSO Factor for Deuteron (S=1) — CRITICAL

**Ptolemy has a defect in the spin-orbit term for S=1 (deuteron, ³He, triton).**

In `WAVELJ`, the L·S coupling constant is:
```fortran
ALS = (J(J+1) - L(L+1) - S(S+1)) / (2 × JSP)
```
where `JSP = 2S`. For S=½ (proton): JSP=1, no change. For S=1 (deuteron): JSP=2, **ALS is half the correct value**.

**Consequence:** Ptolemy's spin-orbit force for deuteron is half-strength.

**Fix for comparison:** When running Ptolemy for (d,d) elastic, multiply your physical VSO by 2:
```
$ Physical VSO = 3.557 MeV (An-Cai 2006)
VSO = 7.114   $ Use 2× in Ptolemy to get correct physics
```

**Does NOT affect unpolarized DCS:** The unpolarized cross section dσ/dΩ is largely insensitive to VSO for S=1 (spin averaging cancels spin-orbit effects). This matters for analyzing powers and spin observables.

### 11.2 Radius Convention Mismatch

- `R0DEFAULT` (Ptolemy native) uses `r0 × (A_t^{1/3} + A_p^{1/3})` for heavy projectiles
- Cleopatra always writes `r0target` → `r0 × A_t^{1/3}`
- **For (d,p): same result** (deuteron A=2 ≤ 2.5 threshold)
- **For (α,α): different by ~14%** at typical nuclear radii

Always specify `r0target` explicitly when matching Cleopatra-generated inputs.

### 11.3 S-Matrix Output Column Confusion

Normal Ptolemy S-matrix output (without `PRINT=1001`) uses the **LX spin-flip basis**, not J=L±½. The columns are NOT directly the JP-split S-matrix elements. Always use `PRINT=1001` to get the physical JP-basis S-matrix.

### 11.4 CM vs Lab Cross Section

The default output table has both CM and Lab cross sections. A common mistake is extracting the Lab column (col 4) for comparison with CM-frame codes. Always use the **C.M. column (col 3)**.

### 11.5 ASYMPTOPIA Auto-Extension

Ptolemy automatically extends the integration grid if the bound state wavefunction doesn't reach its asymptotic form. This triggers warnings like:
```
**** WARNING: SEARCH FOR THE DIF GRIDS WAS STOPPED 7 TIMES BY BOUND STATE ASYMPTOPIA
```
This is normal for loosely bound states. Set `asymptopia=40` or larger to suppress.

---

## 12. Complete Input Examples

### Example 1: Elastic Scattering ⁶⁰Ni(p,p) at 30 MeV

```
HEADER  60Ni(p,p) at 30 MeV  (Koning-Delaroche OMP)
ANGLESTEP=1
ANGLEMIN=1
ANGLEMAX=180

CHANNEL  1H + 60Ni
ELAB=30

V=47.937  R0=1.200  A=0.669
W=2.853   RI0=1.200 AI=0.669
VSI=6.878 RSI0=1.280 ASI=0.550
VSO=5.250 RSO0=1.020 ASO=0.590
VSOI=-0.162 RSOI0=1.020 ASOI=0.590
RC0=1.258

COMPLEXW
ELASTIC SCATTERING
;
END
```

Notes:
- `COMPLEXW` uses volume imaginary (VI = `W`) instead of surface-only
- `VSOI` is negative (imaginary spin-orbit with opposite sign)
- `W` is an alias for `VI`

---

### Example 2: Elastic Scattering ⁶⁰Ni(d,d) at 60 MeV

```
HEADER  60Ni(d,d) at 60 MeV  (An-Cai 2006 OMP)
HEADER  NOTE: VSO=3.557 physical; using VSO=3.557 (Ptolemy half-strength by design)
ANGLESTEP=1
ANGLEMIN=1
ANGLEMAX=180

CHANNEL  2H + 60Ni
ELAB=60.0

V=81.919   R0=1.150  A=0.768
VI=4.836   RI0=1.330 AI=0.464
VSI=8.994  RSI0=1.373 ASI=0.774
VSO=3.557  RSO0=0.972 ASO=1.011
VSOI=0.000 RSOI0=0.972 ASOI=1.011
RC0=1.303

ELASTIC SCATTERING
;
END
```

> ⚠️ Note on VSO: Ptolemy's VSO for deuteron is half the physical value (see §11.1). The An-Cai value VSO=3.557 MeV is designed to be used **as-is** in Ptolemy (which will halve it internally). To get the correct physics in other codes (Raphael, C++), use VSO=3.557/2 = 1.778 MeV, or apply the 1/(2S) correction explicitly.

---

### Example 3: Elastic Scattering ¹⁴⁸Sm(α,α) at 50 MeV

```
HEADER  148Sm(a,a) at 50 MeV
ANGLESTEP=1
ANGLEMIN=1
ANGLEMAX=180

CHANNEL  4He + 148Sm
ELAB=50

V=65.5  R0=1.427  A=0.671
W=29.8  RI0=1.427 AI=0.671
RC0=1.4

ELASTIC SCATTERING
;
END
```

Notes:
- No spin-orbit (S=0 for α particle)
- Simple volume imaginary (`W` = volume imaginary)
- No `r0target` → uses R0DEFAULT = `r0 × (A_t^{1/3} + A_p^{1/3})` for α (A=4 > 2.5)

---

### Example 4: Transfer ¹⁶O(d,p)¹⁷O(g.s.) at 20 MeV — New Style

```
reset
REACTION: 16O(d,p)17O(5/2+ 0.000) ELAB=20.0
PARAMETERSET dpsb labangles r0target lstep=1 lmin=0 lmax=40 maxlextrap=0 asymptopia=30

$ Deuteron internal wavefunction (zero-range approximation via AV18)
PROJECTILE
wavefunction av18 r0=1 a=0.5 l=0
;

$ Target bound state: neutron in 0d5/2 of 17O (n=0, l=2, j=5/2, BE=4.143 MeV)
TARGET
nodes=0 l=2 jp=5/2
v=60.0 r0=1.10 a=0.65
rc0=1.30
writestep=0.1
;

$ Incoming: d + 16O (An-Cai 2006 OMP at Elab=20 MeV)
$ Note: VSO=7.114 = 2 × 3.557 to compensate Ptolemy deuteron SO bug
INCOMING
V=88.9546    R0=1.1489   A=0.7508
VI=2.3480    RI0=1.3446  AI=0.6030
VSI=10.2180  RSI0=1.3943 ASI=0.6872
VSO=7.1140   RSO0=0.9720 ASO=1.0110
VSOI=0.0000  RSOI0=0.0000 ASOI=0.0000
RC0=1.3030
;

$ Outgoing: p + 17O (Koning-Delaroche OMP at ~22 MeV)
OUTGOING
V=49.5434    R0=1.1462   A=0.6753
VI=2.0611    RI0=1.1462  AI=0.6753
VSI=7.6703   RSI0=1.3016 ASI=0.5275
VSO=5.2956   RSO0=0.9338 ASO=0.5900
VSOI=-0.1059 RSOI0=0.9338 ASOI=0.5900
RC0=1.3030
;

anglemin=0 anglemax=180 anglestep=5
;
writens crosssec
end
```

---

### Example 5: S-Matrix Output in JP Basis

```
HEADER  60Ni(p,p) at 30 MeV — JP-basis S-matrix
ANGLEMIN=10
ANGLEMAX=10
ANGLESTEP=10
PRINT=1001

CHANNEL  1H + 60Ni
ELAB=30

V=47.937  R0=1.200  A=0.669
W=2.853   RI0=1.200 AI=0.669
VSI=6.878 RSI0=1.280 ASI=0.550
VSO=5.250 RSO0=1.020 ASO=0.590
VSOI=-0.162 RSOI0=1.020 ASOI=0.590
RC0=1.258

ELASTIC SCATTERING
;
END
```

Output format:
```
ELASTIC S-MATRIX FOR L =  0, JP =  1/2:  0.84697+- 0.52960E-01 I
ELASTIC S-MATRIX FOR L =  1, JP =  1/2:  0.75302+- 0.55983E-01 I
ELASTIC S-MATRIX FOR L =  1, JP =  3/2:  0.88123+- 0.15632E-01 I
...
```

Extraction (Python):
```python
import re
with open('output.out') as f:
    for line in f:
        m = re.search(r'ELASTIC S-MATRIX FOR L =\s*(\d+), JP =\s*(\S+):\s*([\d.E+-]+)\s*\+-?\s*([\d.E+-]+)\s*I', line)
        if m:
            L, JP, Re, Im = m.groups()
            print(f"L={L} JP={JP}: S = {Re} + {Im}i")
```

---

## 13. Quick Reference Card

### Minimal Elastic Input
```
CHANNEL  A + B
ELAB = E
V=...  R0=...  A=...  [other potentials]
ELASTIC SCATTERING
;
END
```

### Minimal Transfer Input (Old-Style)
```
CHANNEL  a + A
ELAB = E
[incoming OMP]
TRANSFER
JBIGA = Jpi_A

CHANNEL  b + B
ELAB = 0.0
[outgoing OMP]
JBIGB = Jpi_B

BOUNDSTATE n l 2j BE
[bound state WS potential]
;
END
```

### Potential Keywords Cheat Sheet
```
V      R0     A       $ real central (depth, radius param, diffuseness)
VI     RI0    AI      $ imaginary volume
VSI    RSI0   ASI     $ imaginary surface (derivative WS)
VSO    RSO0   ASO     $ real spin-orbit
VSOI   RSOI0  ASOI    $ imaginary spin-orbit
RC0                   $ Coulomb radius param
```

### Output Column for Cross Section
```bash
# Extract CM cross section (col 3, 0-based):
grep -E "^\s+[0-9]" output.out | awk '{print $1, $4}'   # theta_CM, sigma_CM
```

### Common Troubleshooting

| Problem | Fix |
|---------|-----|
| `NO SOLUTION FOUND` for bound state | Check V range; try `V=40..100` bracket |
| Norm warning in bound state | Increase `asymptopia` |
| S-matrix not converged | Increase `LMAX` |
| Wrong DCS shape | Check r0target vs R0DEFAULT |
| Deuteron SO mismatch | Use `VSO × 2` for Ptolemy (§11.1) |
| Lab vs CM confusion | Always use col 3 (C.M.) from output |

---

*Manual compiled from Ptolemy source.mor (ANL), validation calculations, and Cleopatra input conventions. Last updated: 2026-03-17.*
