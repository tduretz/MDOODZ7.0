import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl

file = h5py.File('/home/roman/CLionProjects/MDOODZ7.0/cmake-exec/CollisionIra/result6/Output00100.gzip.h5', 'r')
Vx   = file['VxNodes/Vx']
Vz   = file['VzNodes/Vz']
xc   = file['/Model/xc_coord']
zc   = file['/Model/zc_coord']
data = file['/Model/Params']
Nx  = data[3].astype(int)
Nz  = data[4].astype(int)
Ncx = Nx-1
Ncz = Nz-1

Vx  = np.reshape(Vx,  (Nz+1,Nx)).transpose()
Vz  = np.reshape(Vz,  (Nz,Nx+1)).transpose()

# Create 2D centroid mesh
xc2, zc2 = np.meshgrid(xc,zc)

# Average Vx, Vz abd Txy to cell centers
Vxc  = 0.5*( Vx[0:-1, 1:-1 ] + Vx[1:, 1:-1 ] )
Vzc  = 0.5*( Vz[1:-1, 0:-1 ] + Vz[1:-1, 1: ] )

Vx = Vxc * 3.154e9
Vz = Vzc * 3.154e9

# Step for visualisation of quiver arrows
stp = 9


fig, ax = mpl.subplots(nrows=2, ncols=1, figsize=(6, 6))
mpl.set_cmap('jet')
plot1 =  ax[0].pcolor(xc, zc, Vx.transpose())
ax[0].plot()
cbar1 = mpl.colorbar(plot1, ax=ax[0])
cbar1.set_label('Vx [cm/yr]')

plot2 = ax[1].pcolor(xc, zc, Vz.transpose())
ax[1].plot()
cbar2 = mpl.colorbar(plot2, ax=ax[1])
cbar2.set_label('Vz [cm/yr]')


mpl.show()
