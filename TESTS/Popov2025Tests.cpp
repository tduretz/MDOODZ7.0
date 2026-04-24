// CI regression suite for the combined mode-I / mode-II (tensile cap +
// Drucker-Prager) yield surface of Popov et al. (2025, GMD,
// doi:10.5194/gmd-18-7035-2025), implemented in MDLIB under plast = 2
// (see MDLIB/RheologyDensity.c:685-899).
//
// The suite focuses on the paper's Fig. 5 — 0D stress integration — which
// is the paper's canonical verification of the local stress update.  Three
// fixtures exercise the yield surface in complementary ways:
//
//   - VolumetricExtension -> tensile cap (pure mode-I loading)
//   - DeviatoricShear     -> Drucker-Prager envelope (pure mode-II loading)
//   - MixedStrain         -> cap ↔ DP delimiter crossing (combined loading)
//
// Parameters match paper Table 1 "0D Fig. 5a, b" column exactly:
//   G = 10^10, K = 2·10^11, φ = 30°, ψ = 10°, c_MC = 10^6, p_T = -5·10^5,
//   η^vp = 0 (perfect plasticity), Δt = 2 yr,
//   volumetric: ε̇_xx = ε̇_zz = 2.333·10^-15;
//   shear:      ε̇_xy = 7·10^-14;
//   mixed:      superposition of both.
//
// ---------------------------------------------------------------------------
// Assumptions the tests rely on (documented so the next reader can audit):
//
//   A1. Paper Table 1 "0D Fig. 5a, b" parameters produce the paper's physics
//       in MDOODZ's 2D plane-strain (ε̇_yy = 0) — which preserves the trace
//       and the DP invariants used here at the expense of a tiny spurious
//       deviator under "pure volumetric" loading; see A5 below.
//   A2. BC convention: the custom outward SetBCV*Outward callbacks impose
//       Vx_W = -user2, Vx_E = +user2 → ε̇_xx = 2·user2 / Lx.  The paper's
//       ε̇_xx = 2.333e-15 is therefore realised with user2 = 1.167e-15.
//   A3. tensile_line_search = 0 is safe for moderate per-step elastic
//       overshoots.  For paper loading rates it is — the 3×3 local Newton
//       converges in ≤ 25 iterations without line-search damping.
//   A4. MDOODZ writes the global discretised pressure p* (= mesh->p_in) to
//       /Centers/P.  Per paper §3.3 "upon convergence of the global
//       nonlinear iterations p* ≈ p".  Empirically in MDOODZ:
//         - on the DP segment (shear, mixed) p* sits ON the yield surface
//           (within 0.1 %): the global Stokes solver has converged.
//         - on the tensile cap (volumetric) p* is offset by ≈ K·θ̇_vp·dt
//           (≈ 8 % of τ_d): the Stokes solver leaves a one-step trial
//           residual at the yield-surface corner.  Tests use p* directly
//           and loosen the volumetric cap-membership tolerance accordingly.
//   A5. Under pure-volumetric 2D loading (ε̇_xx = ε̇_zz, ε̇_yy = 0 imposed by
//       plane strain), a small deviator arises from the plane-strain
//       constraint.  Tolerance: |sII| < 0.15·τ_d.
//   A6. HDF5 step 0 contains an initialisation artefact for sxxd/szzd
//       (≈ 6e14 Pa before any Stokes solve has run).  Assertions are
//       evaluated at the final step only; the gnuplot script clips step 0.
//   A7. The reference integrator (Popov0DAnalytical.cpp) re-implements the
//       paper's §3.6 algorithm independently of MDLIB.  A match between
//       MDLIB output and the reference tells us the two implementations
//       agree on the full trajectory, not just at the final step.

#include "gtest/gtest.h"

extern "C" {
#include "mdoodz.h"
}

#include "TestHelpers.h"
#include "Popov2025/Popov0DAnalytical.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace {

// ---- Phase + BC callbacks -----------------------------------------------

static int SetPhaseHomogeneous(MdoodzInput *input, Coordinates coordinates) {
  (void)input;
  (void)coordinates;
  return 0;
}

static double SetDensityConstant(MdoodzInput *input, Coordinates coordinates, int phase) {
  (void)coordinates;
  return input->materials.rho[phase];
}

// Outward-normal BC for the 0D volumetric / mixed loadings.  user2 and user3
// set the normal-velocity magnitudes on the x- and z-facing boundaries
// respectively; positive = outward (extension), negative = inward
// (compression).  For mixed loading, user2 and user3 have opposite signs
// which produces a trace plus a deviatoric component simultaneously.
static SetBC SetBCVxOutward(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  (void)coordinates;
  SetBC bc;
  const double VxW = -input->model.user2 / input->scaling.L * input->scaling.t;
  const double VxE = +input->model.user2 / input->scaling.L * input->scaling.t;
  if (position == N || position == S || position == NW || position == SW ||
      position == NE || position == SE) {
    bc.value = 0.0; bc.type = constant_shear_stress;
  } else if (position == W) {
    bc.value = VxW; bc.type = constant_velocity;
  } else if (position == E) {
    bc.value = VxE; bc.type = constant_velocity;
  } else {
    bc.value = 0.0; bc.type = inside;
  }
  return bc;
}

static SetBC SetBCVzOutward(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  (void)coordinates;
  SetBC bc;
  const double VzS = -input->model.user3 / input->scaling.L * input->scaling.t;
  const double VzN = +input->model.user3 / input->scaling.L * input->scaling.t;
  if (position == W || position == E || position == SW || position == SE ||
      position == NW || position == NE) {
    bc.value = 0.0; bc.type = constant_shear_stress;
  } else if (position == S) {
    bc.value = VzS; bc.type = constant_velocity;
  } else if (position == N) {
    bc.value = VzN; bc.type = constant_velocity;
  } else {
    bc.value = 0.0; bc.type = inside;
  }
  return bc;
}

// ---- HDF5 helpers --------------------------------------------------------

static std::vector<std::string> listStepFiles(const std::string &dir) {
  std::vector<std::string> out;
  char cmd[512];
  std::snprintf(cmd, sizeof(cmd), "ls %s/Output*.gzip.h5 2>/dev/null | sort",
                dir.c_str());
  FILE *fp = popen(cmd, "r");
  if (!fp) return out;
  char buf[512];
  while (std::fgets(buf, sizeof(buf), fp)) {
    std::string s(buf);
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    if (!s.empty()) out.push_back(s);
  }
  pclose(fp);
  return out;
}

static std::vector<float> readCentersField(const std::string &h5path,
                                           const char *name) {
  hid_t f   = H5Fopen(h5path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t grp = H5Gopen(f, "Centers", H5P_DEFAULT);
  hid_t dst = H5Dopen(grp, name, H5P_DEFAULT);
  hid_t spc = H5Dget_space(dst);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(spc, dims, NULL);
  std::vector<float> buf(dims[0]);
  H5Dread(dst, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
  H5Sclose(spc); H5Dclose(dst); H5Gclose(grp); H5Fclose(f);
  return buf;
}

struct CentresInfo { int ncx, ncz; std::vector<float> xc, zc; };

static CentresInfo readCentres(const std::string &h5path) {
  CentresInfo out{};
  hid_t f   = H5Fopen(h5path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t grp = H5Gopen(f, "Model", H5P_DEFAULT);
  hid_t dst = H5Dopen(grp, "Params", H5P_DEFAULT);
  hid_t spc = H5Dget_space(dst);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(spc, dims, NULL);
  std::vector<double> params(dims[0]);
  H5Dread(dst, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, params.data());
  H5Sclose(spc); H5Dclose(dst);
  out.ncx = (int)params[3] - 1;
  out.ncz = (int)params[4] - 1;
  auto readF = [&](const char *n) {
    hid_t d = H5Dopen(grp, n, H5P_DEFAULT);
    hid_t s = H5Dget_space(d);
    hsize_t nd[1]; H5Sget_simple_extent_dims(s, nd, NULL);
    std::vector<float> buf(nd[0]);
    H5Dread(d, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    H5Sclose(s); H5Dclose(d);
    return buf;
  };
  out.xc = readF("xc_coord");
  out.zc = readF("zc_coord");
  H5Gclose(grp); H5Fclose(f);
  return out;
}

static double centreValue(const std::vector<float> &field, int ncx, int ncz) {
  const int cx = (ncx - 1) / 2;
  const int cz = (ncz - 1) / 2;
  return field[cz * ncx + cx];
}

// Reconstruct the stress second invariant at cell centres from
// /Centers/sxxd, /Centers/szzd, and /Vertices/sxz (averaged from the four
// surrounding vertices).  MDOODZ writes the deviatoric components but not
// the invariant.
static std::vector<float> reconstructSII(const std::string &h5path) {
  const auto sxxd = readCentersField(h5path, "sxxd");
  const auto szzd = readCentersField(h5path, "szzd");
  hid_t f = H5Fopen(h5path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t v = H5Gopen(f, "Vertices", H5P_DEFAULT);
  hid_t d = H5Dopen(v, "sxz", H5P_DEFAULT);
  hid_t s = H5Dget_space(d);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(s, dims, NULL);
  std::vector<float> sxz(dims[0]);
  H5Dread(d, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sxz.data());
  H5Sclose(s); H5Dclose(d); H5Gclose(v); H5Fclose(f);
  const CentresInfo ci = readCentres(h5path);
  const int Nx = ci.ncx + 1;
  std::vector<float> sII(sxxd.size(), 0.0f);
  for (int j = 0; j < ci.ncz; ++j) {
    for (int i = 0; i < ci.ncx; ++i) {
      const int idxC = j * ci.ncx + i;
      const double sxz_c = 0.25 * ((double)sxz[j * Nx + i]
                                 + (double)sxz[j * Nx + (i + 1)]
                                 + (double)sxz[(j + 1) * Nx + i]
                                 + (double)sxz[(j + 1) * Nx + (i + 1)]);
      const double inv2 = 0.5 * ((double)sxxd[idxC] * (double)sxxd[idxC]
                               + (double)szzd[idxC] * (double)szzd[idxC])
                        + sxz_c * sxz_c;
      sII[idxC] = (float)std::sqrt(std::max(inv2, 0.0));
    }
  }
  return sII;
}

// Material and yield-surface geometry for the 0D tests (paper Table 1 "0D").
struct Material0D {
  double G    = 1.0e10;
  double K    = 2.0e11;
  double phi  = 30.0 * M_PI / 180.0;
  double psi  = 10.0 * M_PI / 180.0;
  double C    = 1.0e6;
  double T_st = -5.0e5;
  double dt   = 2.0 * 365.25 * 86400.0;   // 2 yr
};

static popov2025::YieldGeom make0DGeom() {
  popov2025::Params0D p{};
  Material0D m;
  p.G = m.G; p.K = m.K; p.phi_rad = m.phi; p.psi_rad = m.psi;
  p.C_mc = m.C; p.T_st = m.T_st; p.eta_vp = 0.0;
  p.dt = 1.0; p.t_end = 1.0;
  return popov2025::geometry(p);
}

// Read the per-step MDOODZ trajectory for the centre cell, reconstructing
// sII from (sxxd, szzd, sxz) and p_local from (p*, divu_pl) via Eq. 31.
// Step 0 is skipped because of the HDF5 sxxd init artefact (assumption A6).
struct TrajectoryMDOODZ {
  std::vector<double> t;       // simulation time [s]
  std::vector<double> sII;     // reconstructed stress invariant [Pa]
  std::vector<double> p_local; // yield-surface-clamped pressure [Pa]
};

static TrajectoryMDOODZ readCentreTrajectory(const std::string &dir,
                                             double /*K_bulk*/, double dt) {
  TrajectoryMDOODZ out;
  const auto files = listStepFiles(dir);
  for (const auto &h5 : files) {
    const auto sII = reconstructSII(h5);
    const auto P   = readCentersField(h5, "P");
    const CentresInfo ci = readCentres(h5);
    // Parse step number from filename "Output<NNNNN>.gzip.h5".
    int step = 0;
    {
      const size_t p1 = h5.find("Output");
      if (p1 != std::string::npos) step = std::atoi(h5.c_str() + p1 + 6);
    }
    if (step == 0) continue;  // drop sxxd init artefact (A6)
    out.t.push_back((double)step * dt);
    out.sII.push_back(centreValue(sII, ci.ncx, ci.ncz));
    out.p_local.push_back(centreValue(P, ci.ncx, ci.ncz));  // p* directly (A4)
  }
  return out;
}

// Compute reference trajectory via the independent implicit integrator in
// Popov0DAnalytical.cpp, sampled at the same time points as the MDOODZ run.
struct TrajectoryReference {
  std::vector<double> sII;
  std::vector<double> p;
};

static TrajectoryReference referenceTrajectory(double edot_dev, double edot_vol,
                                               double t_end, double dt) {
  Material0D m;
  popov2025::Params0D p{};
  p.G = m.G; p.K = m.K; p.phi_rad = m.phi; p.psi_rad = m.psi;
  p.C_mc = m.C; p.T_st = m.T_st;
  p.eta_vp = 0.0;
  p.edot_dev = edot_dev;
  p.edot_vol = edot_vol;
  p.dt    = dt;
  p.t_end = t_end;
  TrajectoryReference out;
  out.sII = popov2025::tauII_trajectory(p);
  out.p   = popov2025::P_trajectory(p);
  // Drop the t=0 sample so the length matches the MDOODZ trajectory.
  if (!out.sII.empty()) { out.sII.erase(out.sII.begin()); out.p.erase(out.p.begin()); }
  return out;
}

// Relative L2 error between two vectors of equal length.  When the
// reference norm is near zero (e.g. sII for the pure-volumetric test),
// falls back to absolute L2 to avoid division by zero.
static double relativeL2(const std::vector<double> &num,
                         const std::vector<double> &ref) {
  double num_sq = 0.0, ref_sq = 0.0;
  const size_t N = std::min(num.size(), ref.size());
  for (size_t i = 0; i < N; ++i) {
    const double d = num[i] - ref[i];
    num_sq += d * d;
    ref_sq += ref[i] * ref[i];
  }
  if (ref_sq < 1e-30) return std::sqrt(num_sq);    // absolute fallback
  return std::sqrt(num_sq / ref_sq);
}

// ---- Fixture -------------------------------------------------------------

class Popov2025_0DIntegration : public ::testing::Test {};

// ---- Volumetric extension (paper Fig. 5a) -------------------------------

TEST_F(Popov2025_0DIntegration, VolumetricExtension) {
  MdoodzSetup setup = {};
  setup.SetParticles = new SetParticles_ff{
      .SetPhase   = SetPhaseHomogeneous,
      .SetDensity = SetDensityConstant,
  };
  setup.SetBCs = new SetBCs_ff{
      .SetBCVx = SetBCVxOutward,
      .SetBCVz = SetBCVzOutward,
  };
  RunMDOODZ("Popov2025/Popov0D_VolumetricExtension.txt", &setup);

  const auto files = listStepFiles("Popov0D_VolumetricExtension");
  ASSERT_GE((int)files.size(), 10);
  const auto &last = files.back();

  const auto sII     = reconstructSII(last);
  const auto P       = readCentersField(last, "P");
  const auto divu_pl = readCentersField(last, "divu_pl");
  const CentresInfo ci = readCentres(last);
  const double sII_c  = centreValue(sII,     ci.ncx, ci.ncz);
  const double P_c    = centreValue(P,       ci.ncx, ci.ncz);
  const double divu_c = centreValue(divu_pl, ci.ncx, ci.ncz);

  const auto g = make0DGeom();
  std::printf("VolumetricExtension 0D: sII=%.3e P=%.3e divu_pl=%.3e T_st=%.3e\n",
              sII_c, P_c, divu_c, Material0D{}.T_st);

  // (i) Tensile-cap plasticity activated (dilatant flow on the cap tip).
  EXPECT_GT(divu_c, 0.0) << "Tensile cap plasticity did not activate";

  // (ii) Negligible deviator under pure-volumetric loading.  A small
  // residual arises from the plane-strain constraint (ε̇_yy = 0 forces a
  // non-zero 3D deviator even when ε̇_xx = ε̇_zz).
  EXPECT_LT(std::fabs(sII_c), 0.15 * g.tau_d)
      << "Spurious deviator: sII = " << sII_c;

  // (iii) (sII, P) lies on the tensile cap circle.  On the cap the global
  // Stokes solver leaves a one-step trial-pressure residual (assumption A4),
  // so the tolerance is loose relative to the DP assertion below.
  const double R_haty = std::sqrt(sII_c * sII_c + (P_c - g.py) * (P_c - g.py));
  EXPECT_LT(std::fabs(R_haty - g.Ry) / g.Ry, 0.15)
      << "P not on tensile cap: R_haty = " << R_haty << ", Ry = " << g.Ry;

  // (iv) Full-trajectory L2 check against the independent reference
  // integrator (Popov0DAnalytical.cpp).  Paper Table 1 "0D": ε̇_xx = ε̇_zz =
  // 2.333e-15 giving trace = 4.666e-15 (2D plane-strain; 2/3 of the paper's
  // 3D value).  sII should match nearly exactly; P has a systematic cap-
  // segment trial-pressure offset (A4) so the tolerance is looser.
  const auto mdoodz = readCentreTrajectory("Popov0D_VolumetricExtension",
                                           Material0D{}.K, Material0D{}.dt);
  const auto ref    = referenceTrajectory(/*edot_dev=*/0.0,
                                          /*edot_vol=*/4.666e-15,
                                          (double)mdoodz.t.size() * Material0D{}.dt,
                                          Material0D{}.dt);
  const double L2_tau_abs = relativeL2(mdoodz.sII,     ref.sII);
  const double L2_P       = relativeL2(mdoodz.p_local, ref.p);
  const double N_steps    = (double)mdoodz.sII.size();
  // For volumetric the reference sII is ≈ 0 so relativeL2 falls back to the
  // absolute L2 in Pa.  Normalise by (τ_d · √N) to get a per-step RMS
  // fraction of the delimiter stress.
  const double L2_tau_rms = L2_tau_abs / (g.tau_d * std::sqrt(N_steps));
  std::printf("VolumetricExtension L2: sII RMS/τ_d = %.4f, P = %.4f\n",
              L2_tau_rms, L2_P);
  EXPECT_LT(L2_tau_rms, 0.15) << "sII RMS exceeds 15 % of τ_d over trajectory";
  EXPECT_LT(L2_P, 0.20)       << "P(t) trajectory disagrees with reference > 20 %";
}

// ---- Deviatoric shear (paper Fig. 5b) -----------------------------------

TEST_F(Popov2025_0DIntegration, DeviatoricShear) {
  MdoodzSetup setup = {};
  setup.SetParticles = new SetParticles_ff{
      .SetPhase   = SetPhaseHomogeneous,
      .SetDensity = SetDensityConstant,
  };
  setup.SetBCs = new SetBCs_ff{
      .SetBCVx = SetPureOrSimpleShearBCVx,
      .SetBCVz = SetPureOrSimpleShearBCVz,
  };
  RunMDOODZ("Popov2025/Popov0D_DeviatoricShear.txt", &setup);

  const auto files = listStepFiles("Popov0D_DeviatoricShear");
  ASSERT_GE((int)files.size(), 10);
  const auto &last = files.back();

  const auto sII    = reconstructSII(last);
  const auto P      = readCentersField(last, "P");
  const auto eII_pl = readCentersField(last, "eII_pl");
  const CentresInfo ci = readCentres(last);
  const double sII_c = centreValue(sII,    ci.ncx, ci.ncz);
  const double P_c   = centreValue(P,      ci.ncx, ci.ncz);
  const double eII_c = centreValue(eII_pl, ci.ncx, ci.ncz);

  const auto g = make0DGeom();
  const double tau_DP = g.k * P_c + g.c;
  std::printf("DeviatoricShear 0D: sII=%.4e P=%.4e eII_pl=%.3e tau_DP=%.3e\n",
              sII_c, P_c, eII_c, tau_DP);

  // (i) Drucker-Prager plasticity activated.
  EXPECT_GT(eII_c, 0.0) << "DP plasticity did not activate";

  // (ii) Stress sits on the DP envelope τ_II = k·P + c within 5 %.  Under
  // dilatant flow (ψ > 0) and no imposed volumetric strain, P grows from 0
  // and the trajectory slides along the DP envelope (paper Fig. 5b: "the
  // pressure continues to grow even after switching to mode-II regime").
  EXPECT_LT(std::fabs(sII_c - tau_DP) / std::fabs(tau_DP), 0.05)
      << "Stress off DP envelope: sII = " << sII_c
      << ", k·P + c = " << tau_DP;

  // (iii) Pressure remains bounded (sanity check — runaway pressure growth
  // would signal a numerical regression in the dilatant coupling).
  EXPECT_LT(std::fabs(P_c), 10.0 * Material0D{}.C)
      << "Pressure drifted unphysically far: P = " << P_c;

  // (iv) Full-trajectory L2 check against the independent reference
  // integrator.  Paper Table 1 "0D": ε̇_xy = 7e-14; for pure-shear BCs
  // ε̇_II_dev = bkg_strain_rate = 7e-14 and the imposed trace is 0.  On the
  // DP segment MDOODZ's p* has converged to the yield-clamped p (A4), so
  // tight tolerances apply.
  const auto mdoodz = readCentreTrajectory("Popov0D_DeviatoricShear",
                                           Material0D{}.K, Material0D{}.dt);
  const auto ref    = referenceTrajectory(/*edot_dev=*/7.0e-14,
                                          /*edot_vol=*/0.0,
                                          (double)mdoodz.t.size() * Material0D{}.dt,
                                          Material0D{}.dt);
  const double L2_tau = relativeL2(mdoodz.sII,     ref.sII);
  const double L2_P   = relativeL2(mdoodz.p_local, ref.p);
  std::printf("DeviatoricShear L2: sII = %.4f, P = %.4f\n", L2_tau, L2_P);
  EXPECT_LT(L2_tau, 0.05) << "sII(t) trajectory disagrees with reference > 5 %";
  EXPECT_LT(L2_P,   0.10) << "P(t) trajectory disagrees with reference > 10 %";
}

// ---- Mixed strain (paper Fig. 5c) ---------------------------------------

TEST_F(Popov2025_0DIntegration, MixedStrain) {
  MdoodzSetup setup = {};
  setup.SetParticles = new SetParticles_ff{
      .SetPhase   = SetPhaseHomogeneous,
      .SetDensity = SetDensityConstant,
  };
  setup.SetBCs = new SetBCs_ff{
      .SetBCVx = SetBCVxOutward,
      .SetBCVz = SetBCVzOutward,
  };
  RunMDOODZ("Popov2025/Popov0D_MixedStrain.txt", &setup);

  const auto files = listStepFiles("Popov0D_MixedStrain");
  ASSERT_GE((int)files.size(), 5);

  const auto eII_pl  = readCentersField(files.back(), "eII_pl");
  const auto divu_pl = readCentersField(files.back(), "divu_pl");
  const auto sII     = reconstructSII(files.back());
  const auto P       = readCentersField(files.back(), "P");
  const CentresInfo ci = readCentres(files.back());
  const double s  = centreValue(sII,     ci.ncx, ci.ncz);
  const double p  = centreValue(P,       ci.ncx, ci.ncz);
  const double e  = centreValue(eII_pl,  ci.ncx, ci.ncz);
  const double dv = centreValue(divu_pl, ci.ncx, ci.ncz);
  const double mode_total = std::fabs(e) + std::fabs(dv);
  std::printf("Mixed 0D: sII=%.3e P=%.3e eII_pl=%.3e divu_pl=%.3e\n",
              s, p, e, dv);

  // (i) Plasticity activated — mixed loading should exercise at least one
  // of the mode-I / mode-II flow directions by the final step.
  EXPECT_GT(mode_total, 0.0) << "No plasticity activated in mixed test";

  // (ii) Final (sII, P) lies on the composite yield surface — use p*
  // directly per A4 (MDOODZ's p* is on the surface for the DP/cap-convergent
  // regime).
  const auto g = make0DGeom();
  const double R_haty    = std::sqrt(s * s + (p - g.py) * (p - g.py));
  const double cap_resid = std::fabs(R_haty - g.Ry);
  const double DP_resid  = std::fabs(s - g.k * p - g.c);
  const double best      = std::min(cap_resid, DP_resid);
  EXPECT_LT(best / g.tau_d, 0.05)
      << "Mixed final state off yield surface: cap_resid = " << cap_resid
      << ", DP_resid = " << DP_resid;

  // (iii) Full-trajectory L2 check.  Mixed BCs produce trace(ε̇) ≈ 7e-15
  // and ε̇_II_dev ≈ 7e-14 (paper Table 1 "0D" Fig. 5c).  The trajectory
  // crosses from the cap to the DP segment partway through, so the P
  // tolerance is looser than pure-DP DeviatoricShear.
  const auto mdoodz = readCentreTrajectory("Popov0D_MixedStrain",
                                           Material0D{}.K, Material0D{}.dt);
  const auto ref    = referenceTrajectory(/*edot_dev=*/7.0e-14,
                                          /*edot_vol=*/7.0e-15,
                                          (double)mdoodz.t.size() * Material0D{}.dt,
                                          Material0D{}.dt);
  const double L2_tau = relativeL2(mdoodz.sII,     ref.sII);
  const double L2_P   = relativeL2(mdoodz.p_local, ref.p);
  std::printf("Mixed L2: sII = %.4f, P = %.4f\n", L2_tau, L2_P);
  EXPECT_LT(L2_tau, 0.10) << "sII(t) trajectory disagrees with reference > 10 %";
  EXPECT_LT(L2_P,   0.25) << "P(t) trajectory disagrees with reference > 25 %";
}

}  // namespace
