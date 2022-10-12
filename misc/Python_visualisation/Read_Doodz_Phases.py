import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl

file   = h5py.File('/Users/romankulakov/CLionProjects/MDOODZ70/cmake-exec/BuggySet/Output00000.gzip.h5', 'r')
phases = file['/VizGrid/compo_hr']
xv_ph  = file['/VizGrid/xviz_hr']
zv_ph  = file['/VizGrid/zviz_hr']
data   = file['/Model/Params']

# Achtung the phase map has a twice finer resolution
xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
Ncx_ph = xc_ph.shape[0]
Ncz_ph = zc_ph.shape[0]
phases = np.reshape(phases,(Ncz_ph,Ncx_ph))

mpl.clf()
mpl.pcolor(xc_ph, zc_ph, phases)
mpl.set_cmap('jet')
mpl.colorbar()
mpl.tight_layout()
mpl.show()
mpl.title('T')
