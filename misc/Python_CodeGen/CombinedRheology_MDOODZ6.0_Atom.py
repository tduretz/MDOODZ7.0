from sympy  import *

eta_cst, Eii, eta_ve, eta_el = symbols('eta_cst, Eii, eta_ve, eta_el')
C_pwl, n_pwl = symbols('C_pwl, n_pwl')
C_exp, n_exp, ST = symbols('C_exp, n_exp, ST')
C_lin, n_lin, m_lin, d1 = symbols('C_lin, n_lin, m_lin, d1')
G, gam, pg, lam, cg= symbols('G, gam, pg, lam, cg')
Tii          = 2.0 * eta_ve * Eii
Eii_cst      = Tii/2.0/eta_cst
Eii_pwl      = C_pwl * pow(Tii, n_pwl    );
Eii_exp      = C_exp * pow(Tii, ST+n_exp );
Eii_lin      = C_lin * pow(Tii, n_lin) * pow(d1,-m_lin);
f            = Eii_cst + Eii_pwl + Eii_exp + Eii_lin;


f.diff(eta_ve).subs(d1,'d1').subs(Eii_pwl, 'Eii_pwl').subs(Eii_exp, 'Eii_exp').subs(Eii_lin, 'Eii_lin')


Tii = sqrt(1/2*(Txx**2 + Tzz**2) + Txz**2)
eta_pwl     = (2*C_pwl)**(-1)*Tii**(1-n_pwl)
eta_lin     = (2*C_lin)**(-1)*Tii**(1-n_lin)*d1**m_lin
eta_exp     = (2*C_exp)**(-1)*Tii**(1-(ST+n_exp))
# eta_lin     = (2*C_lin)**(-1))_* pow( Eii, 1.0/n_lin - 1.0 ) * pow(d, m_lin/n_lin);
display(eta_pwl)
eta_ve  = (eta_el**-1 + eta_pwl**-1 + 0*eta_exp**-1 + 0*eta_lin**-1  + 0*eta_cst**-1)**-1
Exx,Ezz,Exz=symbols('Exx,Ezz,Exz')
Txx,Tzz,Txz,P=symbols('Txx,Tzz,Txz,P')


print(eta_ve)

detadTxx = eta_ve.diff(Txx).subs(eta_ve, 'eta_ve').subs(Tii, 'Tii').subs(Tii**2, 'Jii').subs('Jii', 'Tii**2')
detadTzz = eta_ve.diff(Tzz).subs(eta_ve, 'eta_ve').subs(Tii, 'Tii').subs(Tii**2, 'Jii').subs('Jii', 'Tii**2')
detadTxz = eta_ve.diff(Txz).subs(eta_ve, 'eta_ve').subs(Tii, 'Tii').subs(Tii**2, 'Jii').subs('Jii', 'Tii**2')
detadP   = eta_ve.diff(P).subs(eta_ve, 'eta_ve').subs(Tii, 'Tii').subs(Tii**2, 'Jii').subs('Jii', 'Tii**2')

for i in range(1):
   print(detadTxx)
   print(detadTzz)
   print(detadTxz)
   print(detadP)
display(detadTxx)
# d1           = exp(log( G *gam/(lam*(1.0/cg)* Tii *(Eii_pwl + Eii_exp)*pg))/(1.0+pg))

# This is for the progressive reaction of Pips
eta_el  = symbols('eta_el')
eta_pwl = symbols('eta_pwl')
eta_ve  = (eta_el**-1 + eta_pwl**-1 + 0*eta_exp**-1 + 0*eta_lin**-1  + 0*eta_cst**-1)**-1

detadeta_pwl = eta_ve.diff(eta_pwl).subs(eta_ve,'eta_ve')
display(detadeta_pwl)
Xreac  = symbols('Xreac')
pre_factor, T, R = symbols('pre_factor, T, R ')
npwlprod, Eapwlprod, Apwlprod = symbols('npwlprod, Eapwlprod, Apwlprod')
npwlreac, Eapwlreac, Apwlreac = symbols('npwlreac, Eapwlreac, Apwlreac')
a1      = npwlprod + 1.0
a2      = npwlprod + 1.0
f1      = 1-Xreac
f2      = Xreac
sum_down= f1*a1+f2*a2
n_pwl   = (f1*a1*npwlreac  + f2*a2*npwlprod ) / sum_down
Q_pwl   = (f1*a1*Eapwlreac + f2*a2*Eapwlprod) / sum_down
Prod1   = Apwlreac**(f1*a1/sum_down) * Apwlprod**(f2*a2/sum_down);
sum_n   = f1*npwlreac/(npwlreac + 1.0) + f2*npwlprod/(npwlprod + 1.0);
Prod2   = (npwlreac/(npwlreac + 1.0))**(f1*a1*npwlreac/sum_down) * (npwlprod/(npwlprod+1.0))**(f2*a2*npwlprod/sum_down);
A_pwl   = Prod1 * sum_n**(-n_pwl) * Prod2;
F_pwl   = 1.0/6.0*2.0**(1.0/n_pwl) * 3.0**((n_pwl-1.0)/2.0/n_pwl);
B_pwl   = pre_factor * F_pwl * pow(A_pwl,-1.0/n_pwl) * exp( (Q_pwl)/R/n_pwl/T );
C_pwl   = pow(2.0*B_pwl, -n_pwl);
eta_pwl = (2*C_pwl)**(-1)*Tii**(1-n_pwl)
eta_pwl = eta_pwl.simplify()
display(eta_pwl.diff(Xreac))
