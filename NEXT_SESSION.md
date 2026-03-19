# Next Session Briefing — FR-DWBA INELDC Debug
_Generated: 2026-03-19 ~16:50 UTC after Session 4_

## 🎯 Current Status

### What's DONE and VALIDATED ✅
| Component | Status |
|---|---|
| CUBMAP | ✅ validated vs Fortran |
| BSPROD ITYPE=1 | ✅ validated |
| SFROMI 9-J | ✅ validated |
| BETCAL + AMPCAL | ✅ validated (<0.1% DCS) |
| XSECTN + MXZ fix | ✅ validated (inject_smat_test exact match) |
| DCS pipeline end-to-end | ✅ exact Ptolemy CM match when S-matrix injected |

### What's BROKEN ❌
**INELDC — the 6D radial integral.** This is the ONLY remaining bug.

Known errors (pre-9J S-matrix, LSQPOL enabled, vs Ptolemy `print=2`):
- Li=0, Lo=2, Lx=2, JPO=3/2: |S| = −1.64%
- Li=0, Lo=2, Lx=2, JPO=5/2: |S| = +0.12%
- Large-LDEL entries (Li≠Lo) worse: up to 41% for (Li=5,Lo=3,JPO=5/2)

## 📁 Key Files

```
Cpp_AI/
├── fr_o16dp.cpp              # main entry point for full run
├── inject_smat_test.cpp      # S-matrix injection test (VALIDATED ✅)
├── src/dwba/
│   ├── ineldc.cpp            # ← THIS IS WHERE THE BUG IS
│   ├── xsectn.cpp            # validated (has debug prints to clean up)
│   ├── wavelj.cpp            # distorted wave solver
│   └── bound.cpp             # kinematics + bound state
├── Raphael_AI/
│   ├── o16dp_gs_cm1deg.out   # TRUE Ptolemy CM DCS reference (no labangles)
│   └── o16dp_smat_print2.out # Ptolemy pre-9J S-matrix (631 entries) + LAB DCS
```

## 🔨 Build Commands

```bash
cd /home/node/working/ptolemy_2019/Cpp_AI

# Ptolemy-mode (use this for all validation runs):
g++ -O2 -std=c++17 -Iinclude -DHAVE_INELDC_FR fr_o16dp.cpp \
    src/dwba/bound.cpp src/dwba/setup.cpp src/dwba/wavelj.cpp \
    src/dwba/grdset.cpp src/dwba/ineldc.cpp src/dwba/ineldc_zr.cpp \
    src/dwba/a12.cpp src/dwba/xsectn.cpp \
    src/dwba/math_utils.cpp src/dwba/potential_eval.cpp \
    src/dwba/rcwfn.cpp src/elastic/elastic.cpp \
    src/dwba/av18_potential.cpp \
    src/input/Isotope.cpp src/input/Potentials.cpp \
    -DINTERP_PTOLEMY -o fr_o16dp_pto

# inject_smat_test (DCS pipeline test — fast, no INELDC):
g++ -O2 -std=c++17 -Iinclude -DHAVE_INELDC_FR inject_smat_test.cpp \
    src/dwba/bound.cpp src/dwba/setup.cpp src/dwba/wavelj.cpp \
    src/dwba/grdset.cpp src/dwba/ineldc.cpp src/dwba/ineldc_zr.cpp \
    src/dwba/a12.cpp src/dwba/xsectn.cpp \
    src/dwba/math_utils.cpp src/dwba/potential_eval.cpp \
    src/dwba/rcwfn.cpp src/elastic/elastic.cpp \
    src/dwba/av18_potential.cpp \
    src/input/Isotope.cpp src/input/Potentials.cpp \
    -DINTERP_PTOLEMY -o inject_smat_test
```

## 🧹 First Task (cleanup)

Remove debug prints from `xsectn.cpp`:
- `[XSectn] START`
- `[XSectn] SOSWT branch`
- `[SFROMI] elem ...`
- `[XSectn] BETCAL done, BETAS size=...`
- `[XSectn] koffs_map precomputed, size=...`

## 🔬 Main Task: Debug INELDC

### Strategy (in order of priority)

**Option A — Single-element integrand trace:**
For ONE element (Li=0, Lo=2, Lx=2, JPI=1, JPO=3): add prints inside the
`(IU, IV, KPH)` loop to dump integrand value at key (ra, rb, phi_ab) points.
Compare to Ptolemy `print=9` output (if it exists for this reaction).

**Option B — Run Ptolemy with print=9:**
Try adding `print=9` to the Ptolemy input and run to get integrand dumps.
Ptolemy input file: `Raphael_AI/o16dp_gs.in` or similar.

**Option C — 1D phi_ab test:**
Fix ra=rb=4 fm, integrate over phi_ab only. Compare C++ vs Fortran.

**Option D — Convergence study:**
Vary NPDIF (currently 40), NPSUM (currently 40), GAMDIF (currently 5).
Check if the ~1.6% JPO=3/2 error is quadrature resolution or a systematic formula bug.

### Key INELDC parameters (current values in ineldc.cpp)
```
NPDIF = 40    (outer radial quadrature points, ra integration)
NPSUM = 40    (inner radial quadrature points, rb integration)  
NPSUMI = 42   (validated as optimal in Session 3)
GAMDIF = 5    (range parameter for outer quadrature)
Lmax = 40     (partial wave cutoff)
```

### Quantum numbers for 16O(d,p)17O g.s.
```
JA = 2   (deuteron spin × 2)
JB = 1   (proton spin × 2)
JBT = 5  (j_neutron_in_17O × 2 = 2×2.5)
JBP = 1  (j_neutron_in_deuteron × 2 = 2×0.5)
lT = 2   (d-orbital)
lP = 0   (s-orbital in deuteron)
Lx = 2   (only transfer multipolarity)
```

### Ptolemy CM reference (o16dp_gs_cm1deg.out, 1° steps)
```
0°:   35.594 mb/sr
1°:   35.551
5°:   34.552
10°:  31.810
30°:   8.036
45°:   1.420
60°:   2.532
90°:   1.628
120°:  0.916
150°:  0.567
180°:  0.689
```

## 📊 Git State
- Latest commit: `b242ad2` (inject_smat_test validated)
- Branch: master, ~38 commits ahead of origin
- **Ryan: please push from local when convenient**

## ⚡ Quick Sanity Check
```bash
# Verify inject_smat_test still gives 35.594 at 0°:
./inject_smat_test 2>/dev/null | grep "^0.0"
# Expected: 0.0           3.5594e+01
```
