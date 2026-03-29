# Ptolemy Output Format Reference

> Guide to reading the standard output of the Fortran Ptolemy code.
> The output is structured in sequential sections, printed as the calculation progresses.

---

## Table of Contents

1. [Header & Input Echo](#1-header--input-echo)
2. [Bound State Parameters](#2-bound-state-parameters)
3. [Optical Model Parameters](#3-optical-model-parameters)
4. [Reaction Summary](#4-reaction-summary)
5. [Integration Grid Summary](#5-integration-grid-summary)
6. [Transfer S-Matrix Timing](#6-transfer-s-matrix-timing)
7. [Elastic S-Matrix Table](#7-elastic-s-matrix-table)
8. [Transfer S-Matrix Elements](#8-transfer-s-matrix-elements)
9. [Cross Section Table](#9-cross-section-table)
10. [Total Cross Section](#10-total-cross-section)

---

## 1. Header & Input Echo

```
                                                 P T O L E M Y
 April 2007 version     Computer: Intel Pentium, Redhat Linux  hostname          29 Mar 26  023817.420

0INPUT... REACTION: 206Hg(d,p)207Hg(9/2+ 0.000) ELAB= 14.780
0INPUT... PARAMETERSET dpsb r0target
0INPUT... lstep=1 lmin=0 lmax=30 maxlextrap=0 asymptopia=50
```

- Lines beginning with `0INPUT...` echo the input file verbatim
- The `0` column-1 character is a Fortran carriage control (page break / double space)
- Warnings appear as `0**** WARNING:` lines

---

## 2. Bound State Parameters

Printed separately for projectile and target bound states:

```
        PROJECTILE BOUND STATE PARAMETERS
E =   -2.2246 MEV     KAPPA = 0.23161
 PROJECTILE MASS =   1.00 AMU     TARGET MASS =   1.00 AMU     REDUCED MASS =    469.46 MEV/C**2
 L =  0        0 NODES
 PROJECTILE SPIN =  1/2     TARGET SPIN =  1/2
 J PROJECTILE =  1/2
```

| Field | Meaning |
|-------|---------|
| E | Binding energy (negative, MeV) |
| KAPPA | Bound state wave number $\kappa = \sqrt{2\mu|E|}/\hbar$ (fm⁻¹) |
| L | Orbital angular momentum |
| NODES | Number of radial nodes (n-1 where n is the principal quantum number) |
| J PROJECTILE | Total angular momentum j = l ± 1/2 |

For the AV18 deuteron wavefunction, additional info appears:

```
 THE BOUND WAVEFUNCTION WAS COMPUTED BY THE AV18     LINKULE:
          Argonne v18 deuteron wave function,     S-state probability = 0.9422
          spectroscopic amplitude ("SPAM") = 0.97069
```

The target bound state also shows the **depth search** if the potential depth V is not fixed:

```
 FOR V =   73.75,  E =  -3.34, WAVEFUNCTION HAS  2 NODES BUT  1 ARE DESIRED.
 WILL TRY A NEW V:      40.000
```

This means Ptolemy iterates V until the correct binding energy and node count are obtained.

---

## 3. Optical Model Parameters

Printed for both incoming and outgoing channels:

```
        OPTICAL MODEL SCATTERING FOR THE INCOMING CHANNEL
E LAB =   14.780 MEV,    E CM =   14.637 MEV,     K =  1.1817         WAVELENGTH =  5.3170 FM
 PROJECTILE MASS =   2.00 AMU,     TARGET MASS = 206.00 AMU,     REDUCED MASS =   1857.47 MEV/C**2

POTENTIAL         COUPLING CONS.    RADIUS    DIFFUSENESS    RADIUS PARAMETER
 REAL CENTRAL          96.891        6.7977      0.7930            1.1510        POWER =   1.0000
 VOLUME ABSORPTION     2.0230        7.8077      0.2640            1.3220        POWER =   1.0000
 SURFACE ABSORPTION    10.378        8.0321      0.8970            1.3600
 REAL SPIN-ORBIT       3.5570        5.7406      1.0110            0.9720
 PT&SPHERE  COULOMB    4.6503        7.6954      0.0000            1.3030
ASYMPTOPIA =   30.000 FM
 STEP SIZE =     0.125 FM        8.000 STEPS PER "WAVELENGTH"
```

| Column | Meaning |
|--------|---------|
| COUPLING CONS. | Potential depth V (MeV) |
| RADIUS | Physical radius R = r₀ × A^(1/3) (fm) |
| DIFFUSENESS | Diffuseness parameter a (fm) |
| RADIUS PARAMETER | Reduced radius r₀ (fm) |
| POWER | Form factor power (1.0 = standard Woods-Saxon) |

- **K** = center-of-mass wave number (fm⁻¹)
- **WAVELENGTH** = 2π/K (fm)
- **PT&SPHERE COULOMB** shows the Sommerfeld parameter η as COUPLING CONS.
- **ASYMPTOPIA** = maximum integration radius (fm)
- **STEP SIZE** = Numerov integration step h (fm)

---

## 4. Reaction Summary

```
                    SUMMARY OF THE REACTION
                       206Hg(d,p)207Hg(9/2+ 0.000)
                 BIGA   (   A   ,    B   )  BIGB   ;   X
M (AMU)         206.00     2.00     1.00   207.00     1.00
Z                80        1        1       80        0
E* (MEV)          0.0000   0.0000   0.0000   0.0000   0.0000
J                 0/2      2/2      1/2      9/2      1/2
PARITY             +1       +1       +1       +1       +1
```

The notation is `BIGA(A, B)BIGB; X` where:

| Symbol | Meaning | Example |
|--------|---------|---------|
| BIGA | Target nucleus | ²⁰⁶Hg |
| A | Projectile | d (deuteron) |
| B | Ejectile | p (proton) |
| BIGB | Residual nucleus | ²⁰⁷Hg |
| X | Transferred particle | n (neutron) |

The **J** row shows spins as `2J/2` (Ptolemy doubled-integer convention): `2/2` = spin-1, `9/2` = spin-9/2.

Below this:

```
     BOUND STATE PROPERTIES
           PROJECTILE   TARGET
 E            -2.2246   -3.3445          Q =   1.1199
 JP             1/2       9/2
 L               0         4
 NODES           0         1
 SPEC. AMP.    0.9707    1.0000
 SPEC. FACTOR  0.9422    1.0000
```

| Field | Meaning |
|-------|---------|
| E | Binding energy (MeV) |
| Q | Reaction Q-value = E_proj + E_targ (MeV) |
| JP | Total angular momentum of transferred particle in each vertex |
| SPEC. AMP. | Spectroscopic amplitude $\mathcal{S}^{1/2}$ |
| SPEC. FACTOR | Spectroscopic factor $\mathcal{S} = (\text{SPEC. AMP.})^2$ |

---

## 5. Integration Grid Summary

```
             SUMMARY OF THE INTEGRATION GRIDS
         NUM. PTS.   MAP TYPE   GAMMA    MINIMUM   "MID. PT."   MAXIMUM
 (RI+RO)/2:   57          2       2.00      0.00      16.82      43.80
 RI-RO:       40          1      12.00      0.191        8.183
 COS(PHI):    20          2       0.00        4        -1.00000     0.66302      0.2000
```

The 3D integral over $(r_i, r_o, \phi)$ is mapped to coordinates:

| Coordinate | Variable | Gauss Points |
|------------|----------|-------------|
| (RI+RO)/2 | Sum coordinate U | NPSUM (typically 40-57) |
| RI-RO | Difference coordinate V | NPDIF (typically 40) |
| COS(PHI) | Azimuthal angle | NPPHI (typically 20) |

**MAP TYPE:** 1 = linear, 2 = rational-sinh (concentrates points near MID. PT.)

```
THERE ARE  1600 POINTS IN THE (RI, RO) GRID
THERE ARE  32000 POINTS IN THE COMPLETE 3-D GRID
```

---

## 6. Transfer S-Matrix Timing

```
 TIME IN A12                        0.000 SEC
 IN PHI & A12 LOOP                  0.190 SEC
 INTERPOLATING H'S                  0.006 SEC
 COMPUTING SCATTERING WAVES         0.000 SEC
 INTERPOLATING SCATTERING WAVES     0.004 SEC
 REST OF RI & RO LOOP               0.008 SEC
 TOTAL TIME                         0.209 SEC

   202720000 PASSES WERE MADE THROUGH THE INNERMOST LOOP
```

This is the main computational cost — the INELDC radial integral loop.

---

## 7. Elastic S-Matrix Table

```
               INCOMING ELASTIC                          UNITARITY                          OUTGOING ELASTIC
    L   L' LX  MAGNITUDE   PHASE  COULOMB       ELASTIC  REACTION  RESIDUAL        L   L' LX  MAGNITUDE   PHASE  COULOMB

    0    0  0  0.154383   -1.139    3.264       0.02383   0.01404   0.96212        0    0  0  0.488883   -0.884    1.255

    4    4  0  0.267524   -0.499    7.646       0.07160   0.01403   0.91437        4    4  0  0.739757   -0.034    5.016
    4    4  1  0.005735    1.195                                                   4    4  1  0.035019    4.020
    4    4  2  0.000340    8.273
```

**Left block — Incoming elastic S-matrix:**

| Column | Meaning |
|--------|---------|
| L | Orbital angular momentum |
| L' | = L (elastic) |
| LX | Spin-orbit index: 0 = average, 1 = J=L+S, 2 = J=L-S (for spin-1 deuteron: 0=average, 1=J=L+1, 2=J=L-1) |
| MAGNITUDE | \|S_L\| (should be ≤ 1) |
| PHASE | Nuclear phase shift δ_L (radians) |
| COULOMB | Coulomb phase σ_L (radians) |

**Center block — Unitarity check:**

| Column | Meaning |
|--------|---------|
| ELASTIC | \|S\|² — elastic fraction |
| REACTION | Fraction lost to transfer reactions |
| RESIDUAL | 1 - ELASTIC - REACTION (should be small if all channels included) |

**Right block — Outgoing elastic S-matrix** (same format as incoming).

For LX > 0, only the S-matrix is shown (no unitarity).

---

## 8. Transfer S-Matrix Elements

```
   TRANSFER S-MATRIX ELEMENTS

    L   L'  LX    |S|          PHASE(RAD)      |S|          PHASE(RAD)      |S|          PHASE(RAD)
                           JP = 1/2                    JP = 3/2                    JP = 5/2

    0    4   4  0.34961E-04   -0.543    0.10283E-03   -0.465    0.14028E-03   -0.454
    1    3   4  0.71820E-04   -0.655    0.16736E-03   -0.591    0.21295E-03   -0.573
```

Each row gives the transfer S-matrix element for a specific $(L_i, L_o, L_x)$ combination, broken down by JP (conserved total angular momentum of the system, in half-integer notation `2J/2`):

| Column | Meaning |
|--------|---------|
| L | Incoming partial wave $L_i$ |
| L' | Outgoing partial wave $L_o$ |
| LX | Transferred angular momentum $L_x$ |
| \|S\| | Magnitude of transfer S-matrix element |
| PHASE | Phase of transfer S-matrix element (radians) |
| JP | Total angular momentum (multiple JP values per row) |

The selection rule $L_i + L_o + L_x$ = even must be satisfied.

---

## 9. Cross Section Table

This is the main result — the angular distribution:

```
  C.M.  REACTION     REACTION   LOW L  HIGH L   % FROM
  ANGLE  C.M. MB        /RUTH     %/L   % ERROR  L>LMAX
                                                          INCOMING      OUTGOING    LX = 3    LX = 4    LX = 5
                                                          C.M. MB      /RUTH      /RUTH

  10.00  0.43060      0.000001   25.61    0.00     0.0   0.67407E+06  1.004659    0.989086   0.38807E-03   0.42986   0.35682E-03
                                                                                             0.89703E-06
```

**Main line (columns left to right):**

| Column | Meaning |
|--------|---------|
| C.M. ANGLE | Center-of-mass scattering angle θ (degrees) |
| C.M. MB | **Transfer reaction** dσ/dΩ (mb/sr) — the main result |
| /RUTH | Transfer dσ/dΩ divided by Rutherford cross section |
| LOW L %/L | Percentage of DCS from the lowest partial wave |
| HIGH L % ERROR | Estimated truncation error from L > Lmax |
| % FROM L>LMAX | Percentage contribution from partial waves above Lmax |

**Right-side columns (same line):**

| Column | Meaning |
|--------|---------|
| INCOMING C.M. MB | Elastic scattering dσ/dΩ for the incoming channel (mb/sr) |
| /RUTH (first) | Incoming elastic σ/σ_Ruth |
| /RUTH (second) | Outgoing elastic σ/σ_Ruth (= 1.0 far from resonances) |

**Lx breakdown (continuation line):**

The DCS decomposed by transferred angular momentum $L_x$. For 206Hg(d,p)207Hg with l=4:
- LX = 3, 4, 5, 6 columns show individual contributions
- LX = 4 dominates (since l_T=4, l_P=0 → primary L_x = 4)
- Other L_x values arise from spin-orbit coupling in the 9-J recoupling

---

## 10. Total Cross Section

At the end of the angular distribution:

```
TOTAL:     20.329
                                                         1530.6                  981.32
                                                                                            0.19182E-01    20.284       0.25301E-01
```

| Line | Meaning |
|------|---------|
| TOTAL (first number) | Integrated transfer cross section (mb) |
| Second line | Integrated elastic cross sections: incoming, outgoing (mb) |
| Third line | Integrated transfer cross section broken down by $L_x$ |

---

## Column Control Characters

Ptolemy uses Fortran carriage control in column 1:

| Character | Meaning |
|-----------|---------|
| `0` | Double space (new paragraph) |
| `1` | New page (form feed) |
| `+` | Overprint (continuation on same line — used for multi-line headers) |
| ` ` (space) | Normal single space |

These characters are artifacts of Fortran `WRITE` statements and may appear literally when redirecting to a file.

---

*Generated by Dudu 🐱 — from Ptolemy output analysis*
