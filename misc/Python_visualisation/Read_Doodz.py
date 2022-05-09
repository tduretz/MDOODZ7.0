import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl

file = h5py.File('/Users/romankulakov/CLionProjects/MDOODZ7/cmake-exec/PinchSwellGSE/Output00160.gzip.h5', 'r')
P    = file['Centers/d']
xc   = file['/Model/xc_coord']
zc   = file['/Model/zc_coord']
data = file['/Model/Params']
data[0]
data[1]
data[2]
Nx  = data[3].astype(int)
Nz  = data[4].astype(int)
Ncx = Nx-1
Ncz = Nz-1

P = np.reshape(P,(Ncz,Ncx))

mpl.clf()
#mpl.pcolor(x,y,h.transpose())
mpl.pcolor(xc,zc,P)
mpl.set_cmap('jet')
mpl.colorbar()
# mpl.clim(0,40)
mpl.tight_layout()
#mpl.contour(x,y,h.transpose(), [500])
#mpl.quiver(xv[::10,::10],yv[::10,::10],hx_itp[::10,::10], hy_itp[::10,::10])
mpl.show()
mpl.title('T')
