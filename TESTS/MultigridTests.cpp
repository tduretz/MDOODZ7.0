// Unit tests for the GMG hierarchy + transfer operators
// (Sections 2, 3, 5.2 of openspec/changes/add-gmg-stokes-solver/tasks.md).

#include <cmath>
#include <cstring>
#include <gtest/gtest.h>
#include <vector>

extern "C" {
#include "MultigridLevels.h"
}

// ---------------------------------------------------------------------------
// Section 2 - Grid hierarchy
// ---------------------------------------------------------------------------

TEST(GMG_LevelCount, AutoOn201x201YieldsFiveLevels) {
  // Spec scenario: 201x201 -> {201, 101, 51, 26, 13}.
  EXPECT_EQ(MultigridComputeDefaultLevelCount(201, 201, 0), 5);
}

TEST(GMG_LevelCount, AutoOnOddRectangular) {
  // 100x150 should coarsen while min >= 16: 100,50,25 ; next is 13<16 so stop.
  // Levels: {100x150, 50x75, 25x38, 13x19} -> 4 levels.
  EXPECT_EQ(MultigridComputeDefaultLevelCount(100, 150, 0), 4);
}

TEST(GMG_LevelCount, OverrideClampsToAvailable) {
  // Override 2 -> always 2 levels.
  EXPECT_EQ(MultigridComputeDefaultLevelCount(201, 201, 2), 2);
  // Override too large -> clamped to max_possible.
  int lots = MultigridComputeDefaultLevelCount(201, 201, 999);
  EXPECT_GE(lots, 2);
  EXPECT_LE(lots, 10);
}

TEST(GMG_LevelCount, TinyGridsFloorAtTwo) {
  // Even a very small fine grid has to produce >= 2 levels if possible.
  EXPECT_GE(MultigridComputeDefaultLevelCount(10, 10, 0), 2);
  EXPECT_EQ(MultigridComputeDefaultLevelCount(3,  3,  0), 2);
}

TEST(GMG_Hierarchy, AllocateFreeYieldsCorrectSizes) {
  MultigridHierarchy H;
  ASSERT_EQ(MultigridHierarchyAllocate(&H, 201, 201, 1.0, 1.0, 0), 0);
  EXPECT_EQ(H.n_levels, 5);
  const int expected_Nx[] = {201, 101, 51, 26, 13};
  for (int k = 0; k < H.n_levels; ++k) {
    EXPECT_EQ(H.levels[k].Nx, expected_Nx[k]) << "level " << k;
    EXPECT_EQ(H.levels[k].Nz, expected_Nx[k]) << "level " << k;
  }
  // Physical domain length must be invariant across levels.
  const double Lf = (H.levels[0].Nx - 1) * H.levels[0].dx;
  for (int k = 1; k < H.n_levels; ++k) {
    double Lk = (H.levels[k].Nx - 1) * H.levels[k].dx;
    EXPECT_NEAR(Lk, Lf, 1e-12 * Lf);
  }
  MultigridHierarchyFree(&H);
}

// ---------------------------------------------------------------------------
// Section 3 - Transfer operators
// ---------------------------------------------------------------------------

static MultigridHierarchy MakePair(int Nx, int Nz) {
  MultigridHierarchy H;
  // Force 2 levels so we have one fine + one coarse pair to play with.
  int rc = MultigridHierarchyAllocate(&H, Nx, Nz, 1.0 / (Nx - 1), 1.0 / (Nz - 1), 2);
  EXPECT_EQ(rc, 0);
  return H;
}

TEST(GMG_Restrict, PreservesConstantForPressure) {
  MultigridHierarchy H = MakePair(33, 33);
  MultigridLevel *F = &H.levels[0];
  MultigridLevel *C = &H.levels[1];
  for (int k = 0; k < F->Ncx * F->Ncz; ++k) F->P[k] = 3.14;
  RestrictPressure(F, C, F->P, C->P);
  for (int k = 0; k < C->Ncx * C->Ncz; ++k) {
    EXPECT_NEAR(C->P[k], 3.14, 1e-12);
  }
  MultigridHierarchyFree(&H);
}

TEST(GMG_Restrict, PreservesConstantForVelocity) {
  MultigridHierarchy H = MakePair(33, 33);
  MultigridLevel *F = &H.levels[0];
  MultigridLevel *C = &H.levels[1];
  const int nVx_f = F->Nx  * F->Ncz;
  const int nVz_f = F->Ncx * F->Nz;
  for (int k = 0; k < nVx_f; ++k) F->Vx[k] = -1.25;
  for (int k = 0; k < nVz_f; ++k) F->Vz[k] =  0.75;
  RestrictVelocity(F, C, F->Vx, C->Vx, 0);
  RestrictVelocity(F, C, F->Vz, C->Vz, 1);
  // RestrictVelocity zeros the boundary rows (Dirichlet is invariant under
  // restriction — boundary DOFs are held as identity by StokesApplyA and
  // the Vanka block, so any non-zero value there would fossilise into the
  // coarse RHS). Only interior faces are guaranteed to carry the constant.
  for (int J = 0; J < C->Ncz; ++J)
    for (int I = 1; I < C->Nx - 1; ++I)
      EXPECT_NEAR(C->Vx[MdIdx_Vx(I, J, C->Nx)], -1.25, 1e-12);
  for (int J = 1; J < C->Nz - 1; ++J)
    for (int I = 0; I < C->Ncx; ++I)
      EXPECT_NEAR(C->Vz[MdIdx_Vz(I, J, C->Ncx)], 0.75, 1e-12);
  MultigridHierarchyFree(&H);
}

TEST(GMG_Prolong, ExactlyReproducesBilinearForSmoothField) {
  MultigridHierarchy H = MakePair(33, 33);
  MultigridLevel *F = &H.levels[0];
  MultigridLevel *C = &H.levels[1];
  // Set coarse Vx to a linear function: v = a*x + b*z. Bilinear prolongation
  // of a linear field must be exactly reproduced on the fine grid (to the
  // extent allowed by clamping near the boundary).
  for (int j = 0; j < C->Ncz; ++j) {
    for (int i = 0; i < C->Nx; ++i) {
      double x = (double)i;
      double z = (double)j + 0.5;
      C->Vx[MdIdx_Vx(i, j, C->Nx)] = 2.0 * x - 3.0 * z + 7.0;
    }
  }
  for (int k = 0; k < F->Nx * F->Ncz; ++k) F->Vx[k] = 0.0;
  ProlongateVelocityAdd(C, F, C->Vx, F->Vx, 0);
  // The coarse field stores f(x, z) = 2x - 3z + 7 sampled at physical
  // coarse positions (I, J+0.5). A bilinear prolongation of a linear field
  // must be exact, so the fine Vx at physical (if_, jf+0.5) in coarse units
  // equals f(if_*0.5, (jf+0.5)/2). We compare strictly in interior points
  // where the clamped boundary does not alter the result.
  for (int jf = 2; jf < F->Ncz - 2; ++jf) {
    for (int if_ = 2; if_ < F->Nx - 2; ++if_) {
      double x_coarse = if_ * 0.5;                  // vertex-aligned in x
      double z_coarse = (jf + 0.5) / 2.0;           // physical z in coarse cell units
      double expected = 2.0 * x_coarse - 3.0 * z_coarse + 7.0;
      EXPECT_NEAR(F->Vx[MdIdx_Vx(if_, jf, F->Nx)], expected, 1e-10)
        << "at i=" << if_ << ", j=" << jf;
    }
  }
  MultigridHierarchyFree(&H);
}

TEST(GMG_Restrict, FullWeightingRowSumsToOne) {
  // A full-weighting restriction of a constant yields the same constant
  // (already covered above). Check also that the variance-reducing property
  // holds: the Vx restriction of a random vector has smaller L2 per-dof than
  // the injection. This is a weak check but catches sign/shape regressions.
  MultigridHierarchy H = MakePair(33, 33);
  MultigridLevel *F = &H.levels[0];
  MultigridLevel *C = &H.levels[1];
  for (int k = 0; k < F->Nx * F->Ncz; ++k) {
    F->Vx[k] = std::sin(0.3 * k) + std::cos(0.17 * k);
  }
  RestrictVelocity(F, C, F->Vx, C->Vx, 0);
  double s = 0.0;
  for (int k = 0; k < C->Nx * C->Ncz; ++k) s += C->Vx[k] * C->Vx[k];
  double rms_coarse = std::sqrt(s / (C->Nx * C->Ncz));
  double sF = 0.0;
  for (int k = 0; k < F->Nx * F->Ncz; ++k) sF += F->Vx[k] * F->Vx[k];
  double rms_fine = std::sqrt(sF / (F->Nx * F->Ncz));
  // Full-weighting is a low-pass; coarse RMS should be <= fine RMS.
  EXPECT_LE(rms_coarse, rms_fine * 1.01);
  MultigridHierarchyFree(&H);
}

TEST(GMG_ViscRestrict, HarmonicOnLargeContrast) {
  MultigridHierarchy H = MakePair(33, 33);
  MultigridLevel *F = &H.levels[0];
  MultigridLevel *C = &H.levels[1];
  // Half the domain stiff (eta=1e6), half soft (eta=1).
  for (int j = 0; j < F->Ncz; ++j) {
    for (int i = 0; i < F->Ncx; ++i) {
      F->etan[MdIdx_P(i, j, F->Ncx)] = (i < F->Ncx / 2) ? 1.0 : 1e6;
    }
  }
  for (int k = 0; k < F->Nx * F->Nz; ++k) F->etas[k] = 1.0;
  RestrictViscosity(F, C);
  // At the interface column, coarse etan should be well below the arithmetic
  // average (which would be ~5e5). Harmonic mean of {1,1e6,1,1e6} ~= 2.0.
  int mid_I = C->Ncx / 2 - 1;
  for (int J = 0; J < C->Ncz; ++J) {
    double e = C->etan[MdIdx_P(mid_I, J, C->Ncx)];
    EXPECT_LT(e, 1e5) << "interface coarse etan should be closer to harmonic, got " << e;
  }
  MultigridHierarchyFree(&H);
}

// ---------------------------------------------------------------------------
// Section 5.2 - Active-mask restriction monotonicity
// ---------------------------------------------------------------------------

TEST(GMG_ActiveMask, CoarseActiveIfAnyFineActive) {
  MultigridHierarchy H = MakePair(33, 33);
  MultigridLevel *F = &H.levels[0];
  MultigridLevel *C = &H.levels[1];
  // Mark only a single fine cell inactive -> coarse must still be active there
  // iff any of the other 3 fine children are active.
  // First set a 2x2 fine block all inactive -> coarse parent must be inactive.
  std::memset(F->act_P, 1, F->Ncx * F->Ncz);
  F->act_P[MdIdx_P(4, 4, F->Ncx)] = 0;
  F->act_P[MdIdx_P(5, 4, F->Ncx)] = 0;
  F->act_P[MdIdx_P(4, 5, F->Ncx)] = 0;
  F->act_P[MdIdx_P(5, 5, F->Ncx)] = 0;
  RestrictActiveMask(F, C);
  EXPECT_EQ(C->act_P[MdIdx_P(2, 2, C->Ncx)], 0);

  // If only three of four fine children are inactive, coarse remains active.
  F->act_P[MdIdx_P(5, 5, F->Ncx)] = 1;
  RestrictActiveMask(F, C);
  EXPECT_EQ(C->act_P[MdIdx_P(2, 2, C->Ncx)], 1);
  MultigridHierarchyFree(&H);
}

// ---------------------------------------------------------------------------
// Section 5.1 - BuildActiveMask from MDOODZ BC tag arrays
// ---------------------------------------------------------------------------

TEST(GMG_ActiveMask, FromBCTagsRespectsAirRegion) {
  MultigridHierarchy H = MakePair(17, 17);
  MultigridLevel *L = &H.levels[0];
  const int Nx = L->Nx, Nz = L->Nz;
  const int Ncx = L->Ncx, Ncz = L->Ncz;

  // Build BC tag arrays in MDOODZ layout:
  //   BCu: Nx * (Nz + 1) (l in [0, Nz])
  //   BCv: (Nx + 1) * Nz
  //   BCp: Ncx * Ncz
  std::vector<signed char> BCu((size_t)Nx * (Nz + 1), -1);
  std::vector<signed char> BCv((size_t)(Nx + 1) * Nz, -1);
  std::vector<signed char> BCp((size_t)Ncx * Ncz,      -1);

  // Mark top third as air (tag 30) on the pressure cell centres,
  // and the corresponding Vx/Vz faces as well.
  for (int j = (int)(2 * Ncz / 3); j < Ncz; ++j) {
    for (int i = 0; i < Ncx; ++i) BCp[(size_t)i + (size_t)j * Ncx] = 30;
  }
  for (int l = (int)(2 * Ncz / 3) + 1; l <= Nz; ++l) {
    for (int i = 0; i < Nx; ++i) BCu[(size_t)i + (size_t)l * Nx] = 30;
  }
  for (int j = (int)(2 * Ncz / 3) + 1; j < Nz; ++j) {
    for (int k = 0; k < Nx + 1; ++k) BCv[(size_t)k + (size_t)j * (Nx + 1)] = 30;
  }
  // West boundary: type 0 Dirichlet. These should REMAIN active (identity row).
  for (int l = 0; l <= Nz; ++l) BCu[(size_t)0 + (size_t)l * Nx] = 0;

  BuildActiveMask(L, BCu.data(), BCv.data(), BCp.data());

  // Interior bottom third: active.
  EXPECT_EQ(L->act_P [MdIdx_P (Ncx/2, Ncz/4, Ncx)], 1);
  EXPECT_EQ(L->act_Vx[MdIdx_Vx(Nx/2,  Ncz/4, Nx )], 1);
  EXPECT_EQ(L->act_Vz[MdIdx_Vz(Ncx/2, Nz/4,  Ncx)], 1);
  // Top air region: inactive.
  EXPECT_EQ(L->act_P[MdIdx_P(Ncx/2, Ncz - 1, Ncx)], 0);
  // Dirichlet (type 0) boundary: still considered active.
  EXPECT_EQ(L->act_Vx[MdIdx_Vx(0, Ncz/4, Nx)], 1);
  MultigridHierarchyFree(&H);
}
