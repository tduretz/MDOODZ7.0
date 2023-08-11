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

for step in range(file_start, file_end+1, file_step):
    filename    = home + path + 'Output%05d' % step + '.gzip.h5'
    print(filename)
    file   = h5py.File(filename, 'r')
    T     = file['/Centers/T']
    data  = file['/Model/Params']
    xc    = file['/Model/xc_coord']
    zc    = file['/Model/zc_coord']
    Nx    = data[3].astype(int)
    Nz    = data[4].astype(int)
    Ncx   = Nx-1
    Ncz   = Nz-1

    T   = np.reshape( T, (Ncz, Ncx)) - 0*273.15

    mpl.clf()
    mpl.pcolor(xc, zc, T)
    mpl.set_cmap('jet')
    mpl.colorbar()
    mpl.tight_layout()
    mpl.show()
    mpl.title('T')

    file.close()  
