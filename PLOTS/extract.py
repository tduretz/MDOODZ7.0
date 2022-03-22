import matplotlib.pyplot as plt
import numpy as np
import h5py as h5py

fig, axs = plt.subplots(1, 2)
plt.set_cmap('inferno')
plt.tight_layout()

h5_filename = 'ShearTemplateReference.gzip.h5'
file = h5py.File(h5_filename, 'r')
Vx = file['VxNodes/Vx']
Vz = file['VzNodes/Vz']

xc = file['/Model/xc_coord']
zc = file['/Model/zc_coord']
data = file['/Model/Params']
Nx = data[3].astype(int)
Nz = data[4].astype(int)

Vx = np.reshape(Vx, (Nz + 1, Nx)).transpose()
Vz = np.reshape(Vz, (Nz, Nx + 1)).transpose()

# Create 2D centroid mesh
xc2, zc2 = np.meshgrid(xc, zc)

# Average Vx, Vz abd Txy to cell centers
arr1 = Vx[0:-1, 1:-1]
arr2 = Vx[1:, 1:-1]

Vxc = 0.5 * (Vx[0:-1, 1:-1] + Vx[1:, 1:-1])
Vzc = 0.5 * (Vz[1:-1, 0:-1] + Vz[1:-1, 1:])

# Step for visualisation of quiver arrows
stp = 9

axs[0].quiver(xc2[::stp, ::stp], zc2[::stp, ::stp], Vxc[::stp, ::stp].transpose(), Vzc[::stp, ::stp].transpose())
axs[0].set_title(h5_filename)
axs[0].set_aspect('equal')