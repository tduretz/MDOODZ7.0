from sympy import *

Exx = symbols('Exx')
eta1 = Function('eta1')(Exx)
eta2 = Function('eta2')(Exx)
eta3 = Function('eta3')(Exx)
x1,x2,x3 = symbols('x1,x2,x3')

eta = x1*eta1 + x2*eta2 + x3*eta3
eta = 1/(x1/eta1 + x2/eta2 + x3/eta3)
eta = exp(x1*ln(eta1) + x2*ln(eta2) + x3*ln(eta3))
display(eta)
detadexx = eta.diff(Exx)
display(detadexx)
