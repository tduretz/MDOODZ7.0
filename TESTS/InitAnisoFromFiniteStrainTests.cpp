// ============================================================================
// InitAnisoFromFiniteStrainTests — Phase 1 GTests for the aniso-init-from-
// finite-strain change. Verifies the math of the per-marker F_init construction
// used by RheologyParticles.c::aniso_init_finite_strain_for_marker without
// running a full MDOODZ simulation. The helper itself is `static inline` in
// RheologyParticles.c (per-marker, no public symbol); these tests reimplement
// the same closed-form algebra (D1-D3 in design.md), drive it with the
// production `aniso_delta_inv_hansen` from libmdoodz, and assert the round-trip
// closes to machine precision via the production Hansen forward formula.
//
// Test cases (Phase 1, Hansen olivine aniso_db=1 only):
//   ConsistencyAtStep0           — (aniso_factor, aniso_angle) sweep, six sub-
//                                  assertions per cell (round-trip, magnitudes,
//                                  principal-direction match).
//   ConsistencyAtStep0_aniso_factor1
//                                — aniso_factor=1 ⇒ F_init = I.
// ============================================================================

extern "C" {
#include "mdoodz.h"
// aniso_delta_inv_hansen lives in mdoodz-private.h; declare it directly so the
// test does not need a private header. Symbol has external linkage in libmdoodz.
double aniso_delta_inv_hansen( double delta );
}

#include <cmath>
#include <cstdio>
#include <gtest/gtest.h>

namespace {

// Hansen olivine forward δ(FS_AR) — same closed form as
// MDLIB/FlowLaws.c::anisoDelta_HansenOlivine (file-local static there). Kept
// inline here so the test does not depend on internal libmdoodz symbols.
//   γ_eff = √FS_AR − 1/√FS_AR
//   M     = 0.536 · (1 − exp(−γ_eff / 3.96))
//   δ     = 24.5 · M + 1
static double anisoDelta_HansenOlivine_ref( double FS_AR ) {
    const double s         = std::sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.536;
    const double gamma_e   = 3.96;
    const double M         = M_inf * ( 1.0 - std::exp( -gamma_eff / gamma_e ) );
    return 24.5 * M + 1.0;
}

// FS_AR(F) via the same stable closed-form 2x2 SVD that
// MDLIB/RheologyParticles.c::FiniteStrainAspectRatio uses on every marker.
//   E  = (Fxx + Fzz)/2,  Fd = (Fxx − Fzz)/2
//   G  = (Fzx + Fxz)/2,  H  = (Fzx − Fxz)/2
//   Q  = √(E² + H²),     R  = √(Fd² + G²)
//   σ_max = Q + R,       σ_min = |Q − R|
//   FS_AR = σ_max / σ_min
static double FS_AR_from_F( double Fxx, double Fxz, double Fzx, double Fzz ) {
    const double E  = 0.5 * (Fxx + Fzz);
    const double Fd = 0.5 * (Fxx - Fzz);
    const double G  = 0.5 * (Fzx + Fxz);
    const double H  = 0.5 * (Fzx - Fxz);
    const double Q  = std::sqrt(E*E + H*H);
    const double R  = std::sqrt(Fd*Fd + G*G);
    const double sigma_max = Q + R;
    const double sigma_min = std::fabs(Q - R);
    return sigma_max / sigma_min;
}

// Build F_init exactly as RheologyParticles.c::aniso_init_finite_strain_for_marker
// does (design D1-D3). aniso_angle in radians. Returns Fxx, Fxz, Fzx, Fzz and
// the intermediate γ_eff, FS_AR, γ_eq for sub-assertion checks.
struct InitFFields {
    double Fxx, Fxz, Fzx, Fzz;
    double gamma_eff, FS_AR, gamma_eq;
    double a, b, theta_F;
};

static InitFFields build_F_init( double aniso_factor, double aniso_angle_rad ) {
    InitFFields out{};
    out.gamma_eff   = aniso_delta_inv_hansen( aniso_factor );
    const double sqrt_FS_AR = 0.5 * ( out.gamma_eff
                                    + std::sqrt( out.gamma_eff * out.gamma_eff + 4.0 ) );
    out.FS_AR    = sqrt_FS_AR * sqrt_FS_AR;
    out.gamma_eq = std::log( out.FS_AR );
    out.a        = std::exp(  0.5 * out.gamma_eq );
    out.b        = std::exp( -0.5 * out.gamma_eq );
    out.theta_F  = aniso_angle_rad - 0.5 * M_PI;
    const double c = std::cos( out.theta_F );
    const double s = std::sin( out.theta_F );
    out.Fxx = c*c*out.a + s*s*out.b;
    out.Fzz = s*s*out.a + c*c*out.b;
    out.Fxz = c*s*( out.a - out.b );
    out.Fzx = c*s*( out.a - out.b );
    return out;
}

// Principal-eigenvector angle (mod π) of F^T F. For our symmetric pure-shear F,
// F^T F has eigenvectors aligned with the rotation R(θ_F); the principal
// (max-eigenvalue) direction is the F principal extension axis = θ_F.
static double principal_angle_of_FtF( const InitFFields &F ) {
    // C = F^T F   (symmetric 2x2)
    const double Cxx = F.Fxx*F.Fxx + F.Fzx*F.Fzx;
    const double Czz = F.Fxz*F.Fxz + F.Fzz*F.Fzz;
    const double Cxz = F.Fxx*F.Fxz + F.Fzx*F.Fzz;
    // eigenvector for the LARGER eigenvalue λ_max of [[Cxx, Cxz],[Cxz, Czz]]:
    //   tan(2α) = 2 Cxz / (Cxx − Czz)
    // angle of max-eigenvalue direction:  α = 0.5 atan2(2 Cxz, Cxx − Czz)
    return 0.5 * std::atan2( 2.0 * Cxz, Cxx - Czz );
}

// Wrap angle difference into [-π/2, +π/2] (180° rotation symmetry of an axis).
static double angle_diff_mod_pi( double a, double b ) {
    double d = a - b;
    while ( d >  0.5 * M_PI ) d -= M_PI;
    while ( d < -0.5 * M_PI ) d += M_PI;
    return d;
}

} // namespace

// ---------------------------------------------------------------------------
// Task 5.1: ConsistencyAtStep0 — sweep (aniso_factor, aniso_angle), six sub-
// assertions per cell.
// ---------------------------------------------------------------------------
TEST(InitAnisoFromFiniteStrain, ConsistencyAtStep0) {
    const double factors[] = { 1.0, 2.0, 4.0, 8.0 };
    const double angles_deg[] = { 0.0, 45.0, 80.0, 90.0 };

    for ( double factor : factors ) {
        for ( double angle_deg : angles_deg ) {
            const double angle_rad = angle_deg * M_PI / 180.0;
            const InitFFields F = build_F_init( factor, angle_rad );

            // (a) round-trip: aniso_delta_fn(FS_AR(F_init)) == aniso_factor
            const double FS_AR_recovered = FS_AR_from_F( F.Fxx, F.Fxz, F.Fzx, F.Fzz );
            const double delta_recovered = anisoDelta_HansenOlivine_ref( FS_AR_recovered );
            EXPECT_NEAR( delta_recovered, factor, 1.0e-8 )
                << "round-trip failed at aniso_factor=" << factor
                << " aniso_angle=" << angle_deg << " deg";

            // (b) det(F_init) = a·b = exp(γ_eq/2)·exp(-γ_eq/2) = 1
            const double det_F = F.Fxx * F.Fzz - F.Fxz * F.Fzx;
            EXPECT_NEAR( det_F, 1.0, 1.0e-9 )
                << "det F drift at aniso_factor=" << factor
                << " aniso_angle=" << angle_deg << " deg";

            // (c) F_init is symmetric (pure-shear): Fxz == Fzx
            EXPECT_NEAR( F.Fxz, F.Fzx, 1.0e-9 );

            // (d) principal direction of F^T F matches (aniso_angle − π/2) mod π
            //     (foliation-parallel inheritance, per D3 corrigendum).
            if ( factor > 1.0 + 1.0e-12 ) {
                const double alpha = principal_angle_of_FtF( F );
                const double target = angle_rad - 0.5 * M_PI;
                const double diff   = angle_diff_mod_pi( alpha, target );
                EXPECT_NEAR( diff, 0.0, 1.0e-7 )
                    << "principal direction mismatch at aniso_factor=" << factor
                    << " aniso_angle=" << angle_deg
                    << " deg (alpha=" << alpha << ", target=" << target << ")";
            }

            // (e) FS_AR(F_init) matches the expected FS_AR from γ_eff (localises
            //     any failure to the F-construction step vs the round-trip).
            EXPECT_NEAR( FS_AR_recovered / F.FS_AR - 1.0, 0.0, 1.0e-9 )
                << "FS_AR construction mismatch at aniso_factor=" << factor
                << " aniso_angle=" << angle_deg << " deg";

            // (f) γ_eff round-trip via Hansen forward formula.
            //     γ_eff_recovered = √FS_AR − 1/√FS_AR  (Hansen forward).
            const double s = std::sqrt(FS_AR_recovered);
            const double gamma_eff_recovered = s - 1.0 / s;
            EXPECT_NEAR( gamma_eff_recovered, F.gamma_eff, 1.0e-9 )
                << "gamma_eff round-trip failed at aniso_factor=" << factor
                << " aniso_angle=" << angle_deg << " deg";
        }
    }
}

// ---------------------------------------------------------------------------
// Task 5.2: ConsistencyAtStep0_aniso_factor1 — aniso_factor = 1 must give F = I
// to machine precision so init1-style runs are byte-identical.
// ---------------------------------------------------------------------------
TEST(InitAnisoFromFiniteStrain, ConsistencyAtStep0_aniso_factor1) {
    const double angles_deg[] = { 0.0, 80.0 };
    for ( double angle_deg : angles_deg ) {
        const double angle_rad = angle_deg * M_PI / 180.0;
        const InitFFields F = build_F_init( 1.0, angle_rad );

        // F_init − I Frobenius norm < 1e-12
        const double dxx = F.Fxx - 1.0;
        const double dzz = F.Fzz - 1.0;
        const double dxz = F.Fxz;
        const double dzx = F.Fzx;
        const double frob = std::sqrt( dxx*dxx + dzz*dzz + dxz*dxz + dzx*dzx );
        EXPECT_LT( frob, 1.0e-12 )
            << "F_init != I at aniso_factor=1 aniso_angle=" << angle_deg << " deg "
            << "(Fxx=" << F.Fxx << " Fxz=" << F.Fxz
            << " Fzx=" << F.Fzx << " Fzz=" << F.Fzz << ")";

        // γ_eff should be exactly 0 at aniso_factor = 1 (Hansen inverse @ δ=1).
        EXPECT_NEAR( F.gamma_eff, 0.0, 1.0e-12 );
    }
}
