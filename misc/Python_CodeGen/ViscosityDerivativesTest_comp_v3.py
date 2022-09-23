from sympy import *
import numpy as np

Txx, Tyy, Txy,Tzz = symbols('Txx, Tyy, Txy, Tzz ')
Exx,Eyy,Exy,Ezz   = symbols('Exx,Eyy,Exy,Ezz')
C_pwl, n_pwl  = symbols('C_pwl, n_pwl')
eta_0,eta_el  = symbols('eta_0,eta_el')
Tii           = sqrt(1/2*Txx**2 + 1/2*Tyy**2 + 1/2*Tzz**2 + Txy**2);
Jii           = Tii**2
Eii           = sqrt(1/2*Exx**2 + 1/2*Eyy**2 + 1/2*Ezz**2 + Exy**2)

#---------------------- VISCO-ELASTIC VISCOSITY ----------------------#

# 1. Viscosity expression as function of stress invariant
eta_pwl       = (2*C_pwl)**(-1)*Tii**(1-n_pwl)
display(eta_pwl.subs(1/2*Txx**2 + Txy**2 + 1/2*Tyy**2 + 1/2*Tzz**2,'Jii'))
eta_ve  = (eta_el**-1 + eta_pwl**-1)**-1
print('eta_ve = ' + octave_code(eta_ve.subs(eta_pwl,'eta_pwl')) + ';')
# 2. Derivatives of viscosity w.r.t. strain rate components
detadexx = symbols('detadexx')
detadtxx = symbols('detadtxx')
f        = detadexx - ( detadtxx * (2*eta_0 + 2*detadexx*Exx));
detadexx = solve(f,detadexx)[0]
display(detadexx)
print('detadexx = ' + octave_code(detadexx.subs(1/2*Txx**2 + Txy**2 + 1/2*Tyy**2 + 1/2*(Tzz)**2,'Jii')) + ';')
# 3. Derivatives of viscosity w.r.t. stress components
detadtxx = eta_ve.diff(Txx)
detadtyy = eta_ve.diff(Tyy)
detadtxy = eta_ve.diff(Txy)
detadtzz = eta_ve.diff(Tzz)
detadtxx = detadtxx.subs(eta_ve,'eta_ve').subs(1/2*Txx**2 + Txy**2 + 1/2*Tyy**2 + 1/2*Tzz**2,'Jii').simplify()
detadtyy = detadtyy.subs(eta_ve,'eta_ve').subs(1/2*Txx**2 + Txy**2 + 1/2*Tyy**2 + 1/2*Tzz**2,'Jii').simplify()
detadtxy = detadtxy.subs(eta_ve,'eta_ve').subs(1/2*Txx**2 + Txy**2 + 1/2*Tyy**2 + 1/2*Tzz**2,'Jii').simplify()
detadtzz = detadtzz.subs(eta_ve,'eta_ve').subs(1/2*Txx**2 + Txy**2 + 1/2*Tyy**2 + 1/2*Tzz**2,'Jii').simplify()

display(detadtxx)
display(detadtyy)
display(detadtxy)
for i in range(1):
    print('detadTxx = ' + octave_code(detadtxx) + ';')
    print('detadTyy = ' + octave_code(detadtyy) + ';')
    print('detadTxy = ' + octave_code(detadtxy) + ';')
    print('detadTzz = ' + octave_code(detadtzz) + ';')
#---------------------- VISCO-ELASTO-PLASTIC VISCOSITY ----------------------#

# For visco-elasto-plastic viscosity
C,P,phi,eta_Kel,h_phi = symbols('C,P,phi,eta_Kel,h_phi')
lamdot    = Function('lamdot')(Exx,Eyy,Exy,Ezz,P)
Tau_yield = C + P*sin(phi) + eta_Kel*lamdot
eta_vep   = Tau_yield / (2*Eii)

# 1. Visco-elasto-plastic viscosity
display(eta_vep)
# 1. Derivatives w.r.t strain rate components
display(eta_vep.diff(Exx).subs(lamdot,'lamdot'))
display(eta_vep.diff(Exy).subs(lamdot,'lamdot'))

display(eta_vep.diff(P).subs(lamdot,'lamdot'))

#print('deta_vep.dExx = ' + octave_code(eta_vep.diff(Exx)) + ';')

#print('deta_vep.dExy = ' + octave_code(eta_vep.diff(Exy)) + ';')

# 2. Derivatives of lambda dot w.r.t strain rate components
eta_trial = Function('eta_trial')(Exx,Eyy,Exy,Ezz)
F_trial   = Function('F_trial')(Exx,Eyy,Exy,Ezz,P)
lamdot    = F_trial /(eta_trial + eta_Kel + h_phi)
display(lamdot.diff(Exx))
display(lamdot.diff(Exy))
display(lamdot.diff(P))


# 3. Derivatives of F w.r.t strain rate components
Tii     = 2*eta_trial*Eii
F_trial = Tii - P*sin(phi) - C
F_trial.diff(Exx).subs(eta_trial, 'eta_trial').subs(Eii,'Eii')
F_trial.diff(Exy).subs(eta_trial, 'eta_trial').subs(Eii,'Eii')
F_trial.diff(P).subs(eta_trial, 'eta_trial').subs(Eii,'Eii')
