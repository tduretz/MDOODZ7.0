import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl

file = h5py.File('/home/roman/CLionProjects/MDOODZ7.0/cmake-exec/RiftingChenin/result4/Output00200.gzip.h5', 'r')
P    = file['Centers/T']
Vx   = file['VxNodes/Vx']
Vz   = file['VzNodes/Vz']
Txx  = file['Centers/sxxd']
Tzz  = file['Centers/szzd']
Txz  = file['Vertices/sxz']
xc   = file['/Model/xc_coord']
zc   = file['/Model/zc_coord']
data = file['/Model/Params']
data[0]
data[1]
data[2]
Nx  = data[3].astype(int)
Nz  = data[4].astype(int)
Ncx = Nx-1
Ncz = Nz-1

# Reshape arrays
P   = np.reshape(P,   (Ncz,Ncx)).transpose()
Vx  = np.reshape(Vx,  (Nz+1,Nx)).transpose()
Vz  = np.reshape(Vz,  (Nz,Nx+1)).transpose()
Txx = np.reshape(Txx, (Ncz,Ncx)).transpose()
Tzz = np.reshape(Tzz, (Ncz,Ncx)).transpose()
Txz = np.reshape(Txz, (Nz ,Nx )).transpose()

# Create 2D centroid mesh
xc2, zc2 = np.meshgrid(xc,zc) 

# Average Vx, Vz abd Txy to cell centers
Vxc  = 0.5*( Vx[0:-1, 1:-1 ] + Vx[1:, 1:-1 ] )
Vzc  = 0.5*( Vz[1:-1, 0:-1 ] + Vz[1:-1, 1: ] )
Txzc = 0.25*( Txz[0:-1,0:-1] + Txz[0:-1,1:]  + Txz[1:,0:-1]  + Txz[1:,1:] )

# Compute invariant
Tiic = np.sqrt(0.5*(Txx**2 + Tzz**2) + Txzc**2)

# Step for visualisation of quiver arrows
stp = 9

mpl.clf()
mpl.pcolor(xc, zc, P.transpose())
mpl.set_cmap('inferno')
mpl.colorbar()
mpl.tight_layout()
contours = mpl.contour(xc2, zc2, Tiic.transpose(), [2.0], colors='white')
mpl.clabel(contours, inline=True, fontsize=15)
mpl.quiver( xc2[::stp,::stp], zc2[::stp,::stp], Vxc[::stp, ::stp].transpose(), Vzc[::stp, ::stp].transpose())
mpl.title('P map + Tii = 2.0 contour')
mpl.show()

