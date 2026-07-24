(IA generated)
# DocuRheologyParticles

## InitialiseGrainSizeParticles
Initialises the particle grain-size field from the reference grain size of each phase. This is a simple per-marker setup routine used before the rheology loop begins.

## aniso_init_finite_strain_for_marker
Static helper used only for anisotropic finite-strain initialisation on markers. It constructs a deformation gradient consistent with a target anisotropy factor, phase anisotropy angle, and the calibrated inverse anisotropy law, so that the first finite-strain projection is consistent with the prescribed anisotropic fabric.

## InitialiseAnisoDeltaParticles
Initialises the per-marker anisotropy state for phases using ani_fstrain = 3. For anisotropic phases it sets up the marker deformation gradient and relaxed anisotropy variables so the later finite-strain and relaxation updates are internally consistent.

## AccumulatedStrain
Interpolates strain increments from the grid back to particles and accumulates them on marker history fields. This is used to track total strain and the different mechanism-specific strain contributions over time.

## DeformationGradient
Updates the particle deformation gradient tensor from the velocity gradient field. This is the core Lagrangian fabric-tracking step that drives finite-strain anisotropy, because the particle deformation gradient is later used to compute FS_AR.

## FiniteStrainAspectRatio
Computes the finite-strain aspect ratio on particles from the deformation gradient, using a numerically stable 2x2 SVD-based formula. This is the main anisotropy state update: it then projects FS_AR to the grid as FS_AR_n and FS_AR_s, and for ani_fstrain = 3 it also projects the relaxed anisotropy state to aniso_delta_n and aniso_delta_s.

## UpdateMaxPT
Tracks maximum pressure and temperature reached on each particle. This is not anisotropy-specific, but it supports marker history variables used by the model.

## UpdateParticleDensity
Computes density increments on the grid and transfers them back to particles. This updates the particle density history from the current Eulerian density field.

## UpdateParticlePhi
Transfers melt fraction or equivalent phi increments from the grid back to particles. The particle field is updated incrementally and clamped to the physical range [0, 1].

## UpdateParticleX
Transfers the compositional field X from the grid back to particles. The update is incremental, based on the change between the current and previous grid values, with bounds enforcement to keep X in [0, 1].

## UpdateParticleXpips
Updates the reaction-progress or chemical field X using the current grid state and a diffusion-style update. It includes a more elaborate chemical-diffusion workflow and then transfers the resulting increment back to particles.

## UpdateParticleGrainSize
Transfers grain-size updates from the grid to particles. In the current implementation this is the reverse of the usual particle-to-grid workflow for grain size and is used to keep marker grain size consistent with the grid field.

## UpdateParticleEnergy
Updates particle temperature from grid-based thermal increments. Supports both direct incremental transfer and a subgrid-diffusion split, depending on model settings.

## UpdateParticlePressure
Updates particle pressure from grid-based pressure increments. Like the energy update, it can use a subgrid diffusion split when enabled.

## UpdateParticleStress
Updates particle deviatoric stress from the current grid solution. This is important for anisotropy because, when anisotropy and advection are active, it also rotates the director/fabric state on particles so the stored fabric follows the deformation history.

## UpdateParticleDivThermal
Transfers the thermal divergence increment from the grid back to particles. This is used to keep particle thermal-source history consistent with the Eulerian thermal solve.

## UpdateAdvectionMode
Switches the advection flag on or off based on model time. This is a small control routine used during staged simulations.

## UpdateParticlePhase
Switches particle phases after the configured phase-transition time is reached. This updates marker phase IDs directly and is independent of the anisotropy machinery.
