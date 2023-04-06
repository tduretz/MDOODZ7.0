import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
from os.path import expanduser
home = expanduser("~")
 
# Path relative to your home directory
path ="/REPO/MDOODZ7.0/MDLIB/"

# File numbers
file_start = 80
file_step  = 10
file_end   = 80

viz_d = True

for step in range(file_start, file_end+1, file_step):
    filename    = home + path + 'Output%05d' % step + '.gzip.h5'
    print(filename)
    file   = h5py.File(filename, 'r')
    phases = file['/VizGrid/compo_hr']
    xv_ph  = file['/VizGrid/xviz_hr']
    zv_ph  = file['/VizGrid/zviz_hr']
    data   = file['/Model/Params']

    d     = file['/Centers/d']
    xc    = file['/Model/xc_coord']
    zc    = file['/Model/zc_coord']
    Nx    = data[3].astype(int)
    Nz    = data[4].astype(int)
    Ncx   = Nx-1
    Ncz   = Nz-1
    d     = np.reshape( d, (Ncz, Ncx)) 

    # Achtung the phase map has a twice finer resolution
    xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
    zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
    Ncx_ph = xc_ph.shape[0]
    Ncz_ph = zc_ph.shape[0]
    phases = np.reshape(phases,(Ncz_ph,Ncx_ph))

    if viz_d:
        mpl.clf()
        ax = mpl.subplot(1,1,1)
        mpl.pcolor(xc, zc, np.log10(d/1e-6), vmin=1, vmax=3)
        mpl.set_cmap('jet')
        mpl.colorbar()
        mpl.tight_layout()
        ax.set_xlim(-0.4, 0.4)
        ax.set_ylim(-0.175, 0.175)
        mpl.title('d') 
        mpl.show()
    
    else:
        mpl.clf()
        mpl.pcolor(xc_ph, zc_ph, phases)
        mpl.set_cmap('jet')
        mpl.colorbar()
        mpl.tight_layout()
        mpl.title('Phases')
        mpl.show()

    file.close()                             
