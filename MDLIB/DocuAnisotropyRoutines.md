`aniso_delta_inv_hansen`
Analytic inverse for the Hansen olivine anisotropy calibration. Converts a finite-strain anisotropy value delta back to the effective strain argument used by the saturation law, with clamping for invalid or over-saturated inputs.

`AnisoFactorEvolv`
Computes the anisotropy factor used in the viscosity law from finite-strain aspect ratio, mineral calibration, grain-size coupling, and optional relaxation state. Supports the three ani_fstrain modes: raw finite strain, mineral-specific calibration, and relaxed anisotropy.
__Called in UpdateAnisoFactor__

`DeltaRhoProxy`
Internal helper that maps accumulated dislocation-type strain to an effective driving-force proxy for relaxation kinetics. Returns a bounded, monotonic proxy value used by the relaxation timescale model.

`DeltaRelaxationTau`
Computes the anisotropy relaxation timescale from temperature, relaxation length, strain proxy, and mineral-specific kinetics parameters. Returns the scaled relaxation time used by the temperature-dependent anisotropy relaxation model.

`Y2`
Returns the squared second invariant of a deviatoric stress tensor, including the anisotropy-weighted shear contribution.

`I2`
Returns the squared second invariant of a deviatoric strain-rate tensor.

`Solve2x2`
Solves a 2x2 linear system by explicit matrix inversion. Intended as a small local utility for algebraic updates.

`ViscosityConciseAniso`
Main anisotropic rheology kernel. For one phase and one grid point, computes visco-elasto-plastic viscosity, stress components, strain-rate partitioning, overstress, pressure correction, and optional post-processing outputs while accounting for anisotropy, grain-size effects, plasticity, elasticity, and multiple creep mechanisms.
__Called in NonNewtonianViscosityGridAniso__

`UpdateAnisoFactor`
Updates the anisotropy factor on cell centers and vertices by phase-weighted averaging over the local phase composition. Applies the selected averaging scheme and calls the finite-strain anisotropy evolution law when anisotropy is active.
__Called in MAIN__


`NonNewtonianViscosityGridAniso`
Evaluates anisotropic non-Newtonian viscosity and stress fields on the full grid. Interpolates needed fields, calls the anisotropic rheology kernel for each active phase, accumulates phase-weighted properties, and writes the final grid viscosity and stress updates.
__Called in UpdateNonLinearity in StockesRoutines; itself called in MAIN__ no test on non-linearity


`InitialiseDirectorVector`
Initializes particle director vectors from either a per-marker angle or a phase default angle. Optionally adjusts the angle for polar coordinates, then normalizes the resulting direction vector.

`NormalizeDirector`
Normalizes director vectors on cell centers and vertices, skipping free-surface cells. Ensures the anisotropic orientation field remains unit length after interpolation or update.

`AnisotropicDamage`
Applies the anisotropic damage model based on a work-rate threshold and assigns damaged viscosities on cells and vertices. This routine computes a simple isotropic damage proxy and updates viscosity fields accordingly.

`DamagedVolume`
Returns the damaged-volume fraction as a piecewise polynomial of the damage argument. Used by the anisotropic damage model to map the work-rate proxy into a relative weakened volume.


