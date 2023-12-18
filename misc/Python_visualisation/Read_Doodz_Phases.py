import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl

My     = 1e6*3600*24*365.25

file   = h5py.File('/Volumes/T7/d/aniso9/Output03240.gzip.h5', 'r')
# file   = h5py.File('/Users/tduretz/REPO/MDOODZ7.0/MDLIB/Output00000.gzip.h5', 'r')
phases = file['/VizGrid/compo_hr']
xv_ph  = file['/VizGrid/xviz_hr']
zv_ph  = file['/VizGrid/zviz_hr']
data   = file['/Model/Params']
t      = data[0].astype(int)

# Achtung the phase map has a twice finer resolution
xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
Ncx_ph = xc_ph.shape[0]
Ncz_ph = zc_ph.shape[0]
phases = np.reshape(phases,(Ncz_ph, Ncx_ph))

fig, ax = mpl.subplots()
mpl.pcolor(xc_ph, zc_ph, phases)
mpl.tight_layout()
mpl.title('Phases at time = ' + str(t/My) + ' My')
ax.set_aspect('equal', adjustable='box')
mpl.show()
