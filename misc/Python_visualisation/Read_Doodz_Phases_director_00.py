import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
import os

dir_path = '/Users/romankulakov/CLionProjects/cpplots/cmake-build-debug/hdf5/CheninAniso3/'

directory = os.fsencode(dir_path)

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if not filename.endswith(".h5"):
        continue
    file = h5py.File(f'{dir_path}{filename}', 'r')
    nx   = file['Centers/nx']
    nz   = file['Centers/nz']
    xc   = file['/Model/xc_coord']
    zc   = file['/Model/zc_coord']
    data = file['/Model/Params']
    phases = file['/VizGrid/compo_hr']
    xv_ph  = file['/VizGrid/xviz_hr']
    zv_ph  = file['/VizGrid/zviz_hr']
    P    = file['Centers/P']
    xc   = file['/Model/xc_coord']
    zc   = file['/Model/zc_coord']
    Nx  = data[3].astype(int)
    Nz  = data[4].astype(int)
    Ncx = Nx-1
    Ncz = Nz-1

    # Reshape arrays
    nxc  =  np.reshape(nz, (Ncz, Ncx)).transpose()
    # nxc nullify values if they < 0.1
    # nxc[nxc < 0.1] = 0.0
    nzc  = np.reshape(nx, (Ncz, Ncx)).transpose()
    # nzc[nzc > 0.9] = 0.0
    # nxc nullify values if they > 0.9

    # Create 2D centroid mesh
    xc2, zc2 = np.meshgrid(xc,zc)

    xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
    zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
    Ncx_ph = xc_ph.shape[0]
    Ncz_ph = zc_ph.shape[0]
    phases = np.reshape(phases,(Ncz_ph,Ncx_ph))
    P      = np.reshape(P,   (Ncz,Ncx)).transpose()


    # Step for visualisation of quiver arrows
    stp = 9

    mpl.clf()
    mpl.figure()
    mpl.rcParams['figure.figsize'] = 8, 4
    mpl.tight_layout()
    mpl.quiver(xc2[::stp,::stp], zc2[::stp,::stp], nxc[::stp, ::stp].transpose(), nzc[::stp, ::stp].transpose(), headwidth=1.,
               headaxislength=1.,
               headlength=1.,)
    mpl.title('Pressure map and aniso director')
    mpl.savefig(f'{filename}_vel.png', transparent=True)

