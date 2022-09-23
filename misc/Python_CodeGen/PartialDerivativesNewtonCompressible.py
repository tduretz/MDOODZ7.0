from sympy import *

Bpwl, npwl         = symbols('Bpwl, npwl')
Exx, Eyy, Ezz, Exz = symbols('Exx, Eyy, Ezz, Exz')

# Standard invariant including out-of-plane component
I2      = (1/2*(Exx**2 + Eyy**2 + Ezz**2) + Exz**2 )
eta_pwl = Bpwl * sqrt(I2)**(1/npwl-1)
detadexx_v0 = eta_pwl.diff(Exx).subs(I2,'I2')
detadeyy_v0 = eta_pwl.diff(Eyy).subs(I2,'I2') # out-of-plane
display(detadexx_v0)
display(detadeyy_v0)

# Substitute out-of-plane component
Exx, Eyy, Ezz, Exz = symbols('Exx, Eyy, Ezz, Exz')
Eyy     = -(Exx + Ezz)
I2      = (1/2*(Exx**2 + Eyy**2 + Ezz**2) + Exz**2 )
eta_pwl = Bpwl * sqrt(I2)**(1/npwl-1)
display(I2)
detadexx_v1 = eta_pwl.diff(Exx).subs(I2,'I2')
# detadeyy_v1 = eta_pwl.diff(Eyy).subs(I2,'I2') # out-of-plane
display(detadexx_v1)
# display(detadeyy_v1)
