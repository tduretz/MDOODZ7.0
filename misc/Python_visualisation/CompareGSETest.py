import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
import time
from os.path import expanduser
home = expanduser("~")
 
# Path relative to your home directory
path ="/REPO/MDOODZ7.0/MDLIB/"
# path ="/REPO/MDOODZ7.0/VISUAL_TESTS/GSE_REF/"

# File numbers
file_start = 0
file_step  = 10
file_end   = 100

My    = 1e6*3600*24*365.25

for step in range(file_start, file_end+1, file_step):
    filename    = home + path + 'Output%05d' % step + '.gzip.h5'
    print(filename)
    file   = h5py.File(filename, 'r')
    data   = file['/Model/Params']

    d     = file['/Centers/d']
    xc    = file['/Model/xc_coord']
    zc    = file['/Model/zc_coord']
    t     = data[0].astype(int)
    Lx    = data[1].astype(int)
    Lz    = data[2].astype(int)
    Nx    = data[3].astype(int)
    Nz    = data[4].astype(int)
    Ncx   = Nx-1
    Ncz   = Nz-1
    d     = np.reshape( d, (Ncz, Ncx)) 

    mpl.ion()

    fig, ax = mpl.subplots(layout='constrained')
    # mpl.clf()
    mpl.pcolor(xc, zc, np.log10(d), vmin=-5, vmax=-3)
    mpl.set_cmap('jet')
    mpl.colorbar()
    mpl.title('d at time = ' + str(t/My))     
    ax.set_ylim(bottom=-0.25, top=0.25)
    fig.canvas.draw()
    fig.canvas.flush_events()
    # mpl.pause(0.1)
    # time.sleep(0.1)
    
    mpl.show()
    file.close()                             
