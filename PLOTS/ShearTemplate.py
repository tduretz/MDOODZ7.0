import matplotlib.pyplot as plt
import numpy as np
import h5py as h5py


def save_file(h5_filename):
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
    data[0]
    data[1]
    data[2]
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

    fig, axs = plt.subplots(1, 2)

    axs[0, 0].clf()
    axs[0, 0].pcolor(xc, zc, P.transpose())
    axs[0, 0].set_cmap('inferno')
    axs[0, 0].colorbar()
    axs[0, 0].clim(-2.5, 2.5)
    axs[0, 0].tight_layout()
    contours = axs[0, 0].contour(xc2, zc2, Tiic.transpose(), [2.0], colors='white')
    axs[0, 0].clabel(contours, inline=True, fontsize=15)
    axs[0, 0].quiver(xc2[::stp, ::stp], zc2[::stp, ::stp], Vxc[::stp, ::stp].transpose(), Vzc[::stp, ::stp].transpose())
    axs[0, 0].title(f'{h5_filename}: P map + Tii = 2.0 contour')

    axs[0, 1].clf()
    axs[0, 1].pcolor(xc, zc, P.transpose())
    axs[0, 1].set_cmap('inferno')
    axs[0, 1].colorbar()
    axs[0, 1].clim(-2.5, 2.5)
    axs[0, 1].tight_layout()
    contours = axs[0, 1].contour(xc2, zc2, Tiic.transpose(), [2.0], colors='white')
    axs[0, 1].clabel(contours, inline=True, fontsize=15)
    axs[0, 1].quiver(xc2[::stp, ::stp], zc2[::stp, ::stp], Vxc[::stp, ::stp].transpose(), Vzc[::stp, ::stp].transpose())
    axs[0, 1].title(f'{h5_filename}: P map + Tii = 2.0 contour')

    fig.savefig(f'comparison.png', bbox_inches='tight')


save_file('ShearTemplateReference.gzip.h5')
save_file('ShearTemplateResult.gzip.h5')
