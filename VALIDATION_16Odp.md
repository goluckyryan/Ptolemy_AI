# Validation: ¹⁶O(d,p)¹⁷O(g.s.) — Chapter 5.5

**Date:** 2026-03-15
**Reaction:** ¹⁶O(d,p)¹⁷O(g.s. 5/2+) at Elab = 20 MeV (10 MeV/u)
**Transfer:** neutron, l=2, j=5/2, n=0
**Q-value:** 1.9186 MeV
**Binding energies:** Sn(¹⁷O) = 4.1431 MeV, Sn(d) = 2.2246 MeV

## Optical Parameters Used (Handbook Table 5 & 6)

### Incoming: ¹⁶O + d
| Param | Value |
|-------|-------|
| V | 88.955 MeV |
| r₀ | 1.149 fm |
| a | 0.751 fm |
| W_vol | 2.348 MeV |
| r_i | 1.345 fm |
| a_i | 0.603 fm |
| W_surf | 10.218 MeV |
| r_si | 1.397 fm |
| a_si | 0.687 fm |
| V_so | 3.557 MeV |
| r_so | 0.972 fm |
| a_so | 1.011 fm |
| r_c | 1.303 fm |

### Outgoing: ¹⁷O + p
| Param | Value |
|-------|-------|
| V | 49.544 MeV |
| r₀ | 1.146 fm |
| a | 0.675 fm |
| W_vol | 2.061 MeV |
| r_i | 1.146 fm |
| a_i | 0.675 fm |
| W_surf | 7.670 MeV |
| r_si | 1.302 fm |
| a_si | 0.528 fm |
| V_so | 5.296 MeV |
| r_so | 0.934 fm |
| a_so | 0.590 fm |
| V_soi | -0.106 MeV |
| r_c | 1.419 fm |

### Target Bound State (neutron in ¹⁷O)
- n=0, l=2, j=5/2, binding=4.1431 MeV
- r₀=1.25, a=0.65, V_so=6.0, r_so=1.25, a_so=0.65, r_c=1.25
- Ptolemy V_depth = 52.830 MeV
- C++ V_depth = 52.934 MeV (close, ~0.2% difference from mass table differences)

### Projectile Bound State (n in deuteron)
- n=0, l=0, j=1/2, binding=2.2246 MeV
- Ptolemy: AV18 wavefunction, SPAM=0.97069
- C++: r₀=1.0, a=0.5, V_so=0, r_c=1.2

## Results: CM-Frame Differential Cross Sections (mb/sr)

### Ptolemy Finite-Range (reference)
```
θ_CM    dσ/dΩ (mb/sr)
 0°     46.058
 5°     44.926
10°     41.836
15°     36.786
20°     29.315
25°     20.059
30°     11.137
35°      4.803
40°      1.928
45°      1.712
50°      2.634
55°      3.458
60°      3.641
65°      3.251
70°      2.655
75°      2.177
80°      1.951
```
**Total:** 46.588 mb
**Shape:** Peak at 0° (forward-peaked), minimum near ~40° (l=2 pattern), 
secondary maximum ~60°. Classic l=2 stripping angular distribution.

### C++ DWBA Finite-Range (CM frame)
```
θ_CM    dσ/dΩ (mb/sr)
 0°      6.765
 5°      6.373
10°      5.382
15°      4.181
20°      3.069
25°      2.162
30°      1.466
35°      0.953
40°      0.590
45°      0.347
50°      0.194
55°      0.105
60°      0.057
```
**Shape:** Monotonically decreasing — NO minimum or secondary maximum. 
This is WRONG for l=2 transfer.

### Ptolemy Lab-Frame Cross Sections (for completeness)
```
θ_lab   dσ/dΩ_lab (mb/sr)
 0°     53.868
 5°     52.297
10°     48.029
15°     40.959
20°     30.672
25°     18.852
30°      8.889
35°      3.243
40°      1.775
45°      2.551
50°      3.625
55°      4.001
60°      3.624
```

## Zero-Range Status
- **Ptolemy:** Does NOT support zero-range DWBA (FR only).
- **C++ code:** Currently implements FR DWBA only.
- **Handbook Fig 11:** Shows ZR result — would need DWUCK4 or separate ZR code for comparison.
- **D₀² = 1.4–1.55 × 10⁴ MeV²·fm³** (handbook value)

## Comparison Summary

| Angle | Ptolemy FR (mb/sr) | C++ FR (mb/sr) | Ratio C++/Ptol |
|-------|--------------------|----------------|----------------|
| 0°    | 46.058             | 6.765          | 0.147          |
| 10°   | 41.836             | 5.382          | 0.129          |
| 20°   | 29.315             | 3.069          | 0.105          |
| 30°   | 11.137             | 1.466          | 0.132          |
| 40°   |  1.928             | 0.590          | 0.306          |
| 50°   |  2.634             | 0.194          | 0.074          |
| 60°   |  3.641             | 0.057          | 0.016          |

## Diagnosis

1. **Magnitude:** C++ is ~7-15× too small compared to Ptolemy FR.
2. **Shape:** C++ shows monotonic decrease; Ptolemy shows correct l=2 pattern 
   with minimum near 40° and secondary maximum near 60°.
3. **S-matrix elements:** C++ S-matrix elements differ from Ptolemy by factors of 
   1.2-2× in magnitude and have significant phase errors (~0.5-3 radians).
4. **Root causes (likely):**
   - Transfer integral computation (INELDC/radial integration) has systematic errors
   - The 3D finite-range integration grid/method differs from Ptolemy's optimized approach
   - Coupling coefficient (9-J recoupling) may have normalization issues
   - The projectile bound state WS well vs AV18 can cause ~5% difference, not 7×
5. **Bound state potentials** agree well (V=52.83 vs 52.93 MeV).
6. **Kinematics** agree (Ecm=17.75 MeV, k_in=1.235, k_out=0.950).

## Files Generated
- `16O_dp_ptolemy.in` — Ptolemy input (lab angles)
- `16O_dp_cm.in` — Ptolemy input (CM angles)
- `16O_dp_ptolemy.out` — Ptolemy FR output (lab)
- `16O_dp_cm.out` — Ptolemy FR output (CM)
- `o16dp_handbook.in` — C++ DWBA input (handbook parameters)
- `16O_dp_10MeVu.in` — C++ DWBA input (alternative)

## Conclusion
The C++ DWBA code produces qualitatively incorrect results for ¹⁶O(d,p)¹⁷O:
- Wrong angular shape (missing l=2 characteristic pattern)
- Wrong magnitude (factor ~7-15 too small)

The Ptolemy FR calculation gives a physically reasonable result: forward-peaked l=2 
transfer pattern peaking at 0° with minimum near 40° and secondary maximum near 60°.
Total cross section of 46.6 mb is reasonable for this reaction.

The C++ code's transfer integral (radial integration) and/or angular coupling 
coefficients need debugging. The S-matrix element comparison shows systematic 
errors in both magnitude and phase accumulating with L.
