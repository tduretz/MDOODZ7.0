import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
from os.path import expanduser
home = expanduser("~")
 
# Path relative to your home directory
path ="/REPO/MDOODZ7.0/MDLIB/"

# File numbers
file_start = 0
file_step  = 10
file_end   = 100

viz_X = True

for step in range(file_start, file_end+1, file_step):
    filename    = home + path + 'Output%05d' % step + '.gzip.h5'
    print(filename)
    file   = h5py.File(filename, 'r')
    phases = file['/VizGrid/compo_hr']
    xv_ph  = file['/VizGrid/xviz_hr']
    zv_ph  = file['/VizGrid/zviz_hr']
    data   = file['/Model/Params']

    X     = file['/Centers/X']
    xc    = file['/Model/xc_coord']
    zc    = file['/Model/zc_coord']
    Nx    = data[3].astype(int)
    Nz    = data[4].astype(int)
    Ncx   = Nx-1
    Ncz   = Nz-1
    X     = np.reshape( X, (Ncz, Ncx)) 

    # Achtung the phase map has a twice finer resolution
    xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
    zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
    Ncx_ph = xc_ph.shape[0]
    Ncz_ph = zc_ph.shape[0]
    phases = np.reshape(phases,(Ncz_ph,Ncx_ph))

    if viz_X:
        mpl.clf()
        mpl.pcolor(xc, zc, X)
        mpl.set_cmap('jet')
        mpl.colorbar()
        mpl.tight_layout()
        mpl.show()
        mpl.title('X')     
    else:
        mpl.clf()
        mpl.pcolor(xc_ph, zc_ph, phases)
        mpl.set_cmap('jet')
        mpl.colorbar()
        mpl.tight_layout()
        mpl.show()
        mpl.title('Phases')

    file.close()                             
