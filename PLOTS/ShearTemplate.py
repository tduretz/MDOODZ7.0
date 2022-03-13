import matplotlib.pyplot as plt
import numpy as np
import h5py as h5py

fig, axs = plt.subplots(1, 2)
plt.set_cmap('inferno')
plt.tight_layout()

h5_filename = 'ShearTemplateReference.gzip.h5'
file = h5py.File(h5_filename, 'r')
P = file['Centers/P']
Vx = file['VxNodes/Vx']
Vz = file['VzNodes/Vz']
Txx = file['Centers/sxxd']
Tzz = file['Centers/szzd']
Txz = file['Vertices/sxz']
xc = file['/Model/xc_coord']
zc = file['/Model/zc_coord']
data = file['/Model/Params']
Nx = data[3].astype(int)
Nz = data[4].astype(int)
Ncx = Nx - 1
Ncz = Nz - 1

# Reshape arrays
P = np.reshape(P, (Ncz, Ncx)).transpose()
Vx = np.reshape(Vx, (Nz + 1, Nx)).transpose()
Vz = np.reshape(Vz, (Nz, Nx + 1)).transpose()
Txx = np.reshape(Txx, (Ncz, Ncx)).transpose()
Tzz = np.reshape(Tzz, (Ncz, Ncx)).transpose()
Txz = np.reshape(Txz, (Nz, Nx)).transpose()

# Create 2D centroid mesh
xc2, zc2 = np.meshgrid(xc, zc)

# Average Vx, Vz abd Txy to cell centers
Vxc = 0.5 * (Vx[0:-1, 1:-1] + Vx[1:, 1:-1])
Vzc = 0.5 * (Vz[1:-1, 0:-1] + Vz[1:-1, 1:])
Txzc = 0.25 * (Txz[0:-1, 0:-1] + Txz[0:-1, 1:] + Txz[1:, 0:-1] + Txz[1:, 1:])

# Compute invariant
Tiic = np.sqrt(0.5 * (Txx ** 2 + Tzz ** 2) + Txzc ** 2)

# Step for visualisation of quiver arrows
stp = 9

contours = axs[0].contour(xc2, zc2, Tiic.transpose(), [2.0], colors='white')
axs[0].clabel(contours, inline=True, fontsize=15)
im = axs[0].pcolormesh(xc, zc, P.transpose(), shading='auto')
axs[0].quiver(xc2[::stp, ::stp], zc2[::stp, ::stp], Vxc[::stp, ::stp].transpose(), Vzc[::stp, ::stp].transpose())
axs[0].set_title(h5_filename)
axs[0].set_aspect('equal')
plt.colorbar(im, ax=axs[0], location='bottom')

h5_filename = 'ShearTemplateResult.gzip.h5'
file = h5py.File(h5_filename, 'r')
P = file['Centers/P']
Vx = file['VxNodes/Vx']
Vz = file['VzNodes/Vz']
Txx = file['Centers/sxxd']
Tzz = file['Centers/szzd']
Txz = file['Vertices/sxz']
xc = file['/Model/xc_coord']
zc = file['/Model/zc_coord']
data = file['/Model/Params']
Nx = data[3].astype(int)
Nz = data[4].astype(int)
Ncx = Nx - 1
Ncz = Nz - 1

# Reshape arrays
P = np.reshape(P, (Ncz, Ncx)).transpose()
Vx = np.reshape(Vx, (Nz + 1, Nx)).transpose()
Vz = np.reshape(Vz, (Nz, Nx + 1)).transpose()
Txx = np.reshape(Txx, (Ncz, Ncx)).transpose()
Tzz = np.reshape(Tzz, (Ncz, Ncx)).transpose()
Txz = np.reshape(Txz, (Nz, Nx)).transpose()

# Create 2D centroid mesh
xc2, zc2 = np.meshgrid(xc, zc)

# Average Vx, Vz abd Txy to cell centers
Vxc = 0.5 * (Vx[0:-1, 1:-1] + Vx[1:, 1:-1])
Vzc = 0.5 * (Vz[1:-1, 0:-1] + Vz[1:-1, 1:])
Txzc = 0.25 * (Txz[0:-1, 0:-1] + Txz[0:-1, 1:] + Txz[1:, 0:-1] + Txz[1:, 1:])

# Compute invariant
Tiic = np.sqrt(0.5 * (Txx ** 2 + Tzz ** 2) + Txzc ** 2)

# Step for visualisation of quiver arrows
stp = 9

contours = axs[1].contour(xc2, zc2, Tiic.transpose(), [2.0], colors='white')
axs[1].clabel(contours, inline=True, fontsize=15)
im = axs[1].pcolormesh(xc, zc, P.transpose(), shading='auto')
axs[1].quiver(xc2[::stp, ::stp], zc2[::stp, ::stp], Vxc[::stp, ::stp].transpose(), Vzc[::stp, ::stp].transpose())
axs[1].set_title(h5_filename)
axs[1].set_aspect('equal')
plt.colorbar(im, ax=axs[1], location='bottom')

fig.savefig(f'comparison.png', bbox_inches='tight')
