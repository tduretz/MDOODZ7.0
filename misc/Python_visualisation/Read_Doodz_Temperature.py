import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl

file  = h5py.File('/Users/tduretz/REPO/MDOODZ7.0/MDLIB/Output00000.gzip.h5', 'r')
T     = file['/Centers/T']
data  = file['/Model/Params']
xc    = file['/Model/xc_coord']
zc    = file['/Model/zc_coord']
Nx    = data[3].astype(int)
Nz    = data[4].astype(int)
Ncx   = Nx-1
Ncz   = Nz-1

T   = np.reshape( T, (Ncz, Ncx)) - 273.15

print(np.shape(T))
print(np.shape(xc))
print(np.shape(zc))
mpl.clf()
mpl.pcolor(xc, zc, T)
mpl.set_cmap('jet')
mpl.colorbar()
mpl.tight_layout()
mpl.show()
mpl.title('T')
