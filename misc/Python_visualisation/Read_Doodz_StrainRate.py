import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
from os.path import expanduser
home = expanduser("~")
 
# Path relative to your home directory
path ="/REPO/MDOODZ7.0/RUNS/qcoe_x200/"

# File numbers
file_start = 1700
file_step  = 10
file_end   = 1700

viz_eps = True

for step in range(file_start, file_end+1, file_step):
    filename    = home + path + 'Output%05d' % step + '.gzip.h5'
    print(filename)
    file   = h5py.File(filename, 'r')
    phases = file['/VizGrid/compo_hr']
    xv_ph  = file['/VizGrid/xviz_hr']
    zv_ph  = file['/VizGrid/zviz_hr']
    data   = file['/Model/Params']

    xc    = file['/Model/xc_coord']
    zc    = file['/Model/zc_coord']
    Nx    = data[3].astype(int)
    Nz    = data[4].astype(int)
    Ncx   = Nx-1
    Ncz   = Nz-1

    exxd = file['/Centers/exxd']
    ezzd = file['/Centers/ezzd']
    exz = file['/Vertices/exz']
    exxd = np.reshape(exxd,(Ncz,Ncx))
    ezzd = np.reshape(ezzd,(Ncz,Ncx))
    exz = np.reshape(exz,(Nz,Nx))
    exz_c = .25*(exz[0:-1, 0:-1 ] + exz[1:, 0:-1 ] + exz[0:-1, 1: ] + exz[1:, 1: ])
    eII = (.5*exxd**2 + .5*ezzd**2 + exz_c**2)**.5
    log_eII = np.log10(eII)

    # Achtung the phase map has a twice finer resolution
    xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
    zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
    Ncx_ph = xc_ph.shape[0]
    Ncz_ph = zc_ph.shape[0]
    phases = np.reshape(phases,(Ncz_ph,Ncx_ph))

    if viz_eps:
        fig, ax = mpl.subplots()
        mpl.clf()
        mpl.pcolor(xc, zc, log_eII, vmin=-17, vmax=-13)
        mpl.set_cmap('turbo')
        mpl.colorbar()
        mpl.contour(xc_ph, zc_ph, phases, levels = [0, 1], colors=['w','w'])
        mpl.tight_layout()
        mpl.xlim(-0.00258, 0.00258)
        mpl.ylim(-0.00258, 0.00258)
        ax.set_aspect('equal', adjustable='box') 
        mpl.title('log_eII')     
        mpl.show()
    else:
        mpl.clf()
        mpl.pcolor(xc_ph, zc_ph, phases)
        mpl.set_cmap('jet')
        mpl.colorbar()
        mpl.tight_layout()
        mpl.show()
        mpl.title('Phases')

    file.close()                             
