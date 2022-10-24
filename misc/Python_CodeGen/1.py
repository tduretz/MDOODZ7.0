import numpy                    #here we load numpy library that provides a bunch of useful matrix operations akin to MATLAB
from matplotlib import pyplot   #2D plotting library that we will use to plot our results


#setting the grid
space = 6e6
h_litho = 2e6
nx = 100                     #grid points
dx  = space / (nx - 1)          #grid spacing
nt = 20                     #number of timesteps
nu = 0.3                    #viscosity
diff_coef = 1.5
dt = diff_coef*dx**2           #based on von neumaan stability analysis



#innitial condition

u = numpy.ones(nx)
u[:int(h_litho/dx+1)] = 0.040223999999999996      #Square Wave Profile
u[int(h_litho/dx+1):] = -0.019211462686567173      #Square Wave Profile

pyplot.plot(numpy.linspace(0,space,nx), u, label='Initial Solution')

#discritization

un = numpy.ones(nx)

for n in range(nt+1):               #time marching
    un = u.copy()
    for i in range(1, nx-1):        #Space marching
        u[i] = un[i] + nu * dt/dx**2 *(un[i+1] - 2*un[i] + un[i-1]) #Central Differnece Scheme


pyplot.plot(numpy.linspace(0,space,nx), u, label='Numerical Solution')
pyplot.xlabel('Grid Space')
pyplot.ylabel('Velocity')
pyplot.legend()
pyplot.show()