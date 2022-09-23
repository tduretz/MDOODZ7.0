from sympy import *

rho0, beta, Pt, alpha, T, n = symbols('rho0, beta, Pt, alpha, T, n')
rhoc   = rho0*(1+beta*Pt**n)*(1-alpha*T);
rhoc   = rho0*exp(beta*Pt**n)*exp(-alpha*T);
drhodp = rhoc.diff(Pt).subs(rhoc,'rhoc')
drhodp.simplify()

drho2dp2 = rhoc.diff(Pt).diff(Pt).subs(rhoc,'rhoc')
drho2dp2.simplify()

drho2dp2 = drhodp.diff(Pt).subs(drhodp,'drhodp').subs(rhoc,'rhoc')
drho2dp2.simplify()

g,PN,PS,dy  = symbols('g,PN,PS,dy')
rho_cst     = symbols('rho_cst')
rhoN        = Function('rhoN')(PN)
rhoS        = Function('rhoS')(PS)
rhovy       = 1/2*(rhoN+rhoS)
dPdy        = (PN-PS)/dy
fy          = - dPdy + rhovy*g
#fy = - dPdy + rho_cst*g
fy.diff(PN).subs(Derivative(rhoN),'drhodpN')
fy.diff(PS).subs(Derivative(rhoS),'drhodpS')
# for compressibility
PC,PC0,dt= symbols('PC,PC0,dt')
rhoc     = Function('rhoc')(PC)
drhodP   = Function('drhodP')(PC)
beta_eff = 1/rhoc * drhodP
divV     = 0
fp       = beta_eff*(PC-PC0)/dt + divV;
fp.diff(PC).subs(beta_eff,'beta_eff')
