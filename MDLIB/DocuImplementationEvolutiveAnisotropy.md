# Evolutive anisotropy

Idea : As discussed in **Modelling ductile strain localization with evolutive stochastic rheologies, Tommasi et al, 2026**, under certain conditions, forced viscous deformation results in __anisotropic__ strain weakening. 

*Activation flag is ani-fstrain = 4*

### What does it do

Computes a damaged volume in each cell, and from it, determines the isotropic damage (isotopy viscosity drop), anisotropy factor and director vector angle :
- isotropic damage : iso-strain-rate bound
- anisotropy factor : iso-strain-rate bound over iso-stress ration
- angle : direction of maximum shearing

### Usage in txt file

Currently implemented only for powerlaw viscosity and with n=1 or n=3 with elasticity on. 

- `elastic      = 1`
- `pwlv         = 1 `
- `npwl         = 3.0` ou 1.0
- `ani_fstrain  = 4`


It should be possible to combine that with elasticitym plasticity and other creep mecanism but not tested. 

Without advection ok. With advection ok with `constant_dt = 0`.

### Implementation

Hardcoded parameters :
- `gamma0_dam` = viscosity of damaged material
- `WR_th` = work-rate threshold for damage -> around that value damage is localized otherwise no

Computation of damage and anisotropy factor in `AnisotropyRoutines.c: AnisotropicDamage` called from `AnisotropyRoutines.c: UpdateAnisotropyFactor`, fitted function from microscopic data of the paper in `AnisotropyRoutines.c: DamagedVolume`. Damage and anisotropy factor are computed before the mechanical iterations, from the dissipative work-rate (Wdiss) of precedent time-step, advected after projection on the particles (the field has been added on the particles). 

Computation of the anisotropy angle is done in `RheologyParticles.c:UpdateParticleStress` directly on the particles, after the mechanical iteration, using the stress computed at the same time step. So the angle will be used in next time step.




   