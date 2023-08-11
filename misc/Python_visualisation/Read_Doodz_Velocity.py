import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
from os.path import expanduser
home = expanduser("~")

# Path relative to your home directory
path ="/REPO/MDOODZ7.0/MDLIB/"

# File numbers
file_start = 1
file_step  = 10
file_end   = 1

for step in range(file_start, file_end+1, file_step):
    filename    = home + path + 'Output%05d' % step + '.gzip.h5'
    print(filename)
    file   = h5py.File(filename, 'r')
    Vx   = file['VxNodes/Vx']
    Vz   = file['VzNodes/Vz']
    xc   = file['/Model/xc_coord']
    zc   = file['/Model/zc_coord']
    data = file['/Model/Params']
    phases = file['/VizGrid/compo_hr']
    xv_ph  = file['/VizGrid/xviz_hr']
    zv_ph  = file['/VizGrid/zviz_hr']
    data[0]
    data[1]
    data[2]
    Nx  = data[3].astype(int)
    Nz  = data[4].astype(int)
    Ncx = Nx-1
    Ncz = Nz-1

    # Reshape arrays
    Vx  = np.reshape(Vx,  (Nz+1,Nx)).transpose()
    Vz  = np.reshape(Vz,  (Nz,Nx+1)).transpose()

    # Create 2D centroid mesh
    xc2, zc2 = np.meshgrid(xc,zc)

    # Average Vx, Vz to cell centers
    Vxc  = 0.5*( Vx[0:-1, 1:-1 ] + Vx[1:, 1:-1 ] )
    Vzc  = 0.5*( Vz[1:-1, 0:-1 ] + Vz[1:-1, 1: ] )

    xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
    zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
    Ncx_ph = xc_ph.shape[0]
    Ncz_ph = zc_ph.shape[0]
    phases = np.reshape(phases,(Ncz_ph,Ncx_ph))

    # Step for visualisation of quiver arrows
    stp = 9

    mpl.clf()
    mpl.figure()
    # mpl.pcolor(xc_ph, zc_ph, phases)
    mpl.pcolor(xc, zc, Vxc.transpose())
    mpl.set_cmap('jet')
    mpl.colorbar()
    mpl.tight_layout()
    mpl.quiver(xc2[::stp,::stp], zc2[::stp,::stp], Vxc[::stp, ::stp].transpose(), Vzc[::stp, ::stp].transpose(), color='1')
    mpl.title('Phase map')
    # mpl.savefig(f'{filename}_vel.png')
    mpl.show()

