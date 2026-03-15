// ineldc_zr.cpp — DWBA::InelDcZR() — Zero-Range DWBA transfer integral
//
// In the ZR approximation: V_np(rp)*phi_P(rp) = D0 * delta^3(rp)
// This collapses the 6D finite-range integral to a 1D radial integral:
//
//   I_ZR(Li,Lo,Lx,JPI,JPO) = D0 * A12(phi_T=0, phi_ab=0)
//       * integral_0^inf chi_a(ra) * phi_T(ra) * chi_b*(zr_scale*ra) dra
//
// The result is stored in TransferSMatrix with the same format as InelDc,
// so XSectn() can process it identically.

#include "dwba.h"
#include "math_utils.h"
#include "potential_eval.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

// Constants (same as ineldc.cpp)
static const double HBARC_ZR = 197.32697;
static const double AMU_ZR   = 931.494;

// Stub for InelDc — only compiled when FR ineldc.cpp is not linked
#ifndef HAVE_INELDC_FR
void DWBA::InelDc() {
  std::cerr << "ERROR: InelDc() called in ZR binary — use CalculateZR() instead.\n";
  std::abort();
}
#endif

void DWBA::InelDcZR() {
  // ===========================================================
  // Zero-Range DWBA transfer integral
  // ===========================================================

  // ---- Masses ----
  double mA = Incoming.Target.Mass;
  double ma = Incoming.Projectile.Mass;
  double mb = Outgoing.Projectile.Mass;
  double mB = Outgoing.Target.Mass;

  const double AMU_MEV = 931.494061;

  // True mass of transferred particle
  double mx_kinematic = ma - mb;
  double mx = mx_kinematic + ProjectileBS.BindingEnergy / AMU_MEV;

  // ---- ZR coordinate transformation ----
  // Ptolemy GRDSET convention for stripping:
  //   BRATMS(1) = mx/mb  (neutron/proton ≈ 1.001)
  //   BRATMS(2) = mx/mA  (neutron/target)
  double BRATMS1 = mx / mb;
  double BRATMS2 = mx / mA;
  double denom_c = BRATMS1 + BRATMS2 * (1.0 + BRATMS1);
  double S1_c = (1.0 + BRATMS1) * (1.0 + BRATMS2) / denom_c;
  double T1_c = -(1.0 + BRATMS2) / denom_c;
  double S2_c = (1.0 + BRATMS1) / denom_c;
  // ZR scale: at phi_ab=0, rp=0 → rb = (S2_c/S1_c)*ra
  double zr_scale = S2_c / S1_c;   // ≈ 1/(1+BRATMS2) ≈ A/(A+1)
  // Jacobian from GRDSET
  double JACOB_grdset = S1_c * S1_c * S1_c;

  // ---- ZR normalization constant ----
  const double D0_ZR = -120.1;  // MeV·fm^(3/2)

  // ---- Build bound state channels (same as InelDc) ----

  // Target Bound State: neutron x orbiting target core A
  Channel TgtBS_ch;
  TgtBS_ch.Pot    = TargetBS.Pot;
  TgtBS_ch.Target = Incoming.Target;
  TgtBS_ch.Projectile.Z    = Incoming.Projectile.Z - Outgoing.Projectile.Z;
  TgtBS_ch.Projectile.A    = Incoming.Projectile.A - Outgoing.Projectile.A;
  TgtBS_ch.Projectile.Mass = mx;
  TgtBS_ch.mu   = mA * mx / (mA + mx) / AMU_MEV;
  TgtBS_ch.StepSize = Incoming.StepSize;
  TgtBS_ch.MaxR     = Incoming.MaxR;
  TgtBS_ch.NSteps   = Incoming.NSteps;
  TgtBS_ch.RGrid    = Incoming.RGrid;
  TgtBS_ch.WaveFunction.resize(Incoming.NSteps);
  TgtBS_ch.V_real.resize(Incoming.NSteps);
  TgtBS_ch.V_imag.resize(Incoming.NSteps);
  TgtBS_ch.V_so_real.resize(Incoming.NSteps);
  TgtBS_ch.V_so_imag.resize(Incoming.NSteps);
  TgtBS_ch.V_coulomb.resize(Incoming.NSteps);

  WavSet(TgtBS_ch);
  CalculateBoundState(TgtBS_ch, TargetBS.n, TargetBS.l, TargetBS.j, TargetBS.BindingEnergy);

  // Projectile Bound State: NOT NEEDED for ZR (replaced by D0*delta)
  // But we still need ProjectileBS parameters for ATERM computation.

  // ---- Determine allowed Lx ----
  int lT = TargetBS.l;
  int lP = ProjectileBS.l;
  int LxMin = std::abs(lT - lP);
  int LxMax_bs = lT + lP;

  // ---- Partial wave limits ----
  int Lmax = 15;

  TransferSMatrix.clear();

  // Spin of projectile/ejectile (doubled integers) for J-split DW
  const int JA_dw = 2;  // 2 * spin_deuteron
  const int JB_dw = 1;  // 2 * spin_proton

  // ---- SFROMI kinematic factor ----
  double ka_sf = Incoming.k;
  double kb_sf = Outgoing.k;
  double Ea_sf = Incoming.Ecm;
  double Eb_sf = Outgoing.Ecm;
  double FACTOR_sfromi = 2.0 * std::sqrt(ka_sf * kb_sf / (Ea_sf * Eb_sf));

  double h_zr = Incoming.StepSize;
  int N_zr = TgtBS_ch.NSteps;

  int total_entries = 0;

  for (int Li = 0; Li <= Lmax; ++Li) {
    // J-split incoming distorted waves for all JPI
    int JPI_min = std::abs(2*Li - JA_dw);
    int JPI_max = 2*Li + JA_dw;
    std::map<int, std::vector<std::complex<double>>> chi_a_byJPI;
    for (int JPI = JPI_min; JPI <= JPI_max; JPI += 2) {
      WavElj(Incoming, Li, JPI);
      chi_a_byJPI[JPI] = Incoming.WaveFunction;
    }

    for (int Lo = 0; Lo <= Lmax; Lo++) {
      // J-split outgoing waves
      int JPO_min = std::abs(2*Lo - JB_dw);
      int JPO_max = 2*Lo + JB_dw;
      std::map<int, std::vector<std::complex<double>>> chi_b_byJPO;
      for (int JPO = JPO_min; JPO <= JPO_max; JPO += 2) {
        if (JPO < 1) continue;
        WavElj(Outgoing, Lo, JPO);
        chi_b_byJPO[JPO] = Outgoing.WaveFunction;
      }

      for (int Lx = LxMin; Lx <= LxMax_bs; Lx += 2) {
        // Triangle rule: |Li-Lo| <= Lx <= Li+Lo
        if (Lx < std::abs(Li - Lo)) continue;
        if (Lx > Li + Lo) continue;
        // Parity: (Li+Lo+Lx) must be even
        if ((Li + Lo + Lx) % 2 != 0) continue;

        // Precompute A12 angular coupling for (Li, Lo, Lx)
        std::vector<std::tuple<int,int,double>> A12_terms =
            ComputeA12Terms(Li, Lo, Lx, lT, lP);

        // At phi_T=0, phi_ab=0 (ZR geometry): A12 = sum of all coefficients
        double A12_at_zero = EvalA12(A12_terms, 0.0, 0.0);

        // Skip if angular coupling vanishes
        if (std::abs(A12_at_zero) < 1e-15) continue;

        // J-split loops
        for (auto &[JPI, chi_a] : chi_a_byJPI) {
          for (auto &[JPO, chi_b] : chi_b_byJPO) {

            // Interpolation of chi_b at arbitrary r
            auto interp_chi_b = [&](double r) -> std::complex<double> {
              if (r <= 0 || r >= Outgoing.MaxR) return {0.0, 0.0};
              double idx_f = r / h_zr;
              int ii = static_cast<int>(idx_f);
              if (ii >= (int)chi_b.size() - 1) return {0.0, 0.0};
              double frac = idx_f - ii;
              return chi_b[ii] * (1.0 - frac) + chi_b[ii+1] * frac;
            };

            // ---- 1D ZR radial integral ----
            // I_1D = sum_i chi_a(ra_i) * phi_T(ra_i) * conj(chi_b(zr_scale*ra_i)) * h
            //
            // Wave function storage:
            //   chi_a[i] = u_a(r_i) = r * chi_a(r)  (radial wave u=r*chi)
            //   TgtBS_ch.WaveFunction[i] = phi_T(r_i) = u_T(r)/r  (bound state, already divided by r)
            //   chi_b[i] = u_b(r_i) = r * chi_b(r)
            //
            // The DWBA integral in terms of radial functions is:
            //   I = integral chi_a(ra) * phi_T(ra) * chi_b*(rb) * ra * rb * J_zr * dra
            //
            // In terms of u functions: chi = u/r, so chi_a = u_a/ra, chi_b = u_b/rb
            //   I = integral (u_a/ra) * phi_T(ra) * (u_b*/rb) * ra * rb * J_zr * dra
            //     = integral u_a * phi_T * u_b* * J_zr * dra
            //
            // But we need to be more careful. The 3D delta function collapse gives:
            //   I_3D = integral d^3 r_a chi_a^* * V_delta * chi_b * phi_T
            //        = integral_0^inf integral dOmega chi_a(ra) * phi_T(rx=ra) 
            //                         * chi_b*(rb=zr_scale*ra) * D0 * delta_Jacobian * ra^2 dra
            //
            // For the partial wave decomposition, the angular integral produces the A12 kernel.
            // The radial integral is then:
            //   I_radial = integral chi_a_L(ra) * phi_T_l(ra) * chi_b_L*(zr_scale*ra) * ra^2 dra
            //
            // With u = r*chi stored in WaveFunction arrays:
            //   chi_a_L(ra) = u_a[i] / ra
            //   phi_T(ra) = TgtBS[i]  (already = u_T/r)
            //   chi_b_L(rb) = u_b(rb) / rb = interp(rb) / rb
            //
            // So the radial integral becomes:
            //   I = integral (u_a/ra) * phi_T * (u_b*/rb) * ra^2 * dra
            //     = integral u_a * phi_T * (u_b*/rb) * ra * dra
            //     = integral u_a * phi_T * u_b* * (ra/rb) * dra
            //     = integral u_a * phi_T * u_b* * (1/zr_scale) * dra

            std::complex<double> I_1D(0.0, 0.0);

            for (int i = 1; i < N_zr - 1; ++i) {
              double ra = i * h_zr;
              double rb = zr_scale * ra;

              // u_a(ra) — incoming distorted wave u = r*chi
              std::complex<double> u_a_r = chi_a[i];

              // phi_T(ra) — target bound state (already = u_T/r from CalculateBoundState)
              double phi_T_r = TgtBS_ch.WaveFunction[i].real();

              // u_b*(rb) — outgoing distorted wave at rb (interpolated)
              std::complex<double> u_b_conj = std::conj(interp_chi_b(rb));

              // The 1D ZR radial integral:
              // chi_a * phi_T * chi_b* = (u_a/ra) * phi_T * (u_b*/rb)
              // Multiply by ra^2 dra for the radial measure:
              // = u_a * phi_T * u_b* * ra / rb * dra
              // = u_a * phi_T * u_b* / zr_scale * dra
              I_1D += u_a_r * phi_T_r * u_b_conj / zr_scale * h_zr;
            }

            std::complex<double> Integral = D0_ZR * A12_at_zero * I_1D;

            // Phase factor (same as InelDc): i^(Li+Lo+2*Lx+1)
            int ITEST = ((Li + Lo + 2*Lx + 1) % 4 + 4) % 4;
            std::complex<double> phase_factor;
            switch (ITEST) {
              case 0: phase_factor = { 1.0,  0.0}; break;
              case 1: phase_factor = { 0.0,  1.0}; break;
              case 2: phase_factor = {-1.0,  0.0}; break;
              case 3: phase_factor = { 0.0, -1.0}; break;
            }
            Integral *= phase_factor;

            // SFROMI ATERM computation
            // Ptolemy RACAH(JA,JB,JC,JD,JE,JF) = (-1)^((JA+JB+JC+JD)/2) * {JA/2,JB/2,JE/2; JD/2,JC/2,JF/2}
            // Called as RACAH(2*lBT, 2*jBT, 2*lBP, 2*jBP, 2*jX, 2*Lx)
            // So SixJ = {lBT, jBT, jX; jBP, lBP, Lx}
            // Phase = (-1)^(lBT + jBT + lBP + jBP) [note: jBT, jBP half-integer, sum is integer]
            double ATERM_val = 0.0;
            if (Lx >= std::abs(TargetBS.l - ProjectileBS.l) &&
                Lx <= TargetBS.l + ProjectileBS.l) {
              // Only non-zero when triangle(lT, lP, Lx) holds
              // For lP=0: Lx must equal lT
              if (Lx == TargetBS.l || ProjectileBS.l != 0) {
                double jT_bs = TargetBS.j;      // j of neutron in target (e.g. 1.5 or 2.5)
                double jP_bs = ProjectileBS.j;   // j of neutron in projectile (0.5)
                double jx = 0.5;                 // neutron spin

                // SixJ{lBT, jBT, jX; jBP, lBP, Lx}
                double sj = SixJ((double)TargetBS.l, jT_bs, jx,
                                 jP_bs, (double)ProjectileBS.l, (double)Lx);
                // Phase: (-1)^((2*lBT + 2*jBT + 2*lBP + 2*jBP)/2)
                // = (-1)^(lBT + jBT + lBP + jBP)
                // Using doubled integers: phase_sum = 2*lT + 2*jT + 2*lP + 2*jP (all integers)
                int phase_sum_halved = TargetBS.l + (int)std::round(jT_bs)
                                     + ProjectileBS.l + (int)std::round(jP_bs);
                // For half-integer j: 2*j is odd. (2*lT + 2*jT + 2*lP + 2*jP)/2
                // = lT + jT + lP + jP. Since jT and jP are half-integers, use 2*j to track parity:
                int twoj_sum = 2*TargetBS.l + (int)(2*jT_bs + 0.5) + 2*ProjectileBS.l + (int)(2*jP_bs + 0.5);
                // twoj_sum/2 gives the exponent (always integer since sum of two odds + two evens)
                double sign_val = ((twoj_sum / 2) % 2 == 0) ? 1.0 : -1.0;
                double RACAH_val = sign_val * sj;

                // JBIGA/JBIGB: nuclear spins from DWBA object (generic)
                int JBIGA = (int)std::round(2.0 * SpinTarget);    // 2*J(target)
                int JBIGB = (int)std::round(2.0 * SpinResidual);  // 2*J(residual)

                double TEMP_aterm = std::sqrt((JBIGB + 1.0) / (JBIGA + 1.0));

                // Spectroscopic amplitudes
                double SPAMP = ProjectileWFLoaded ? ProjectileWFSpam : 0.97069;
                double SPAMT = 1.0;

                ATERM_val = TEMP_aterm * std::sqrt(2.0*Lx + 1.0) * SPAMP * SPAMT * RACAH_val;

                // Sign: Ptolemy ATERM sign convention
                // ITEST = JX - JBP + 2*(LBP + LBT) where all doubled
                int JX_doubled = 1;
                int JBP_doubled = (int)(2 * ProjectileBS.j);
                int ITEST_aterm = JX_doubled - JBP_doubled + 2*(ProjectileBS.l + TargetBS.l);
                if ((ITEST_aterm / 2 + 1) % 2 != 0) ATERM_val = -ATERM_val;
              }
            }

            // Total SFROMI factor: FACTOR * |ATERM| / sqrt(2*Li+1)
            double sfromi_norm = FACTOR_sfromi * std::abs(ATERM_val) / std::sqrt(2.0 * Li + 1.0);
            if (ATERM_val < 0) {
              Integral = -Integral;
            }

            // Store in TransferSMatrix
            auto S_elem = Integral * sfromi_norm;

            if (std::abs(S_elem) > 1e-15) {
              TransferSMatrix.push_back({Lx, Li, Lo, JPI, JPO, S_elem});
              total_entries++;
            }

          } // end JPO loop
        } // end JPI loop
      } // end Lx loop
    } // end Lo loop
  } // end Li loop

  std::printf("[ZR] Transfer integrals computed: %d entries.\n", total_entries);
}
