import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
from os.path import expanduser
home = expanduser("~")

My     = 1e6*3600*24*365.25
 
# Path relative to your home directory
path ="/REPO/MDOODZ7.0/MDLIB/"
path ="/REPO/MDOODZ7.0/RUNS/RiverTom/"

# File numbers
file_start = 228
file_step  = 10
file_end   = 228

for step in range(file_start, file_end+1, file_step):
    filename    = home + path + 'Output%05d' % step + '.gzip.h5'
    print(filename)
    file   = h5py.File(filename, 'r')
    
    data   = file['/Model/Params']
    t      = data[0].astype(int)
    xc     = file['/Model/xc_coord']
    zc     = file['/Model/zc_coord']
    Nx     = data[3].astype(int)
    Nz     = data[4].astype(int)
    Ncx    = Nx-1
    Ncz    = Nz-1

    eta_n  = file['/Centers/eta_n']
    eta_n  = np.reshape( eta_n, (Ncz, Ncx))
    rho_n  = file['/Centers/rho_n']
    rho_n  = np.reshape( rho_n, (Ncz, Ncx))
    P      = file['/Centers/P']
    P      = np.reshape( P, (Ncz, Ncx))

    fig, ax = mpl.subplots()
    mpl.clf()
    # mpl.pcolor(xc, zc, np.log10(eta_n), vmin=19, vmax=25)
    # mpl.pcolor(xc, zc, P)
    mpl.pcolor(xc, zc, rho_n)
    mpl.set_cmap('jet')
    mpl.colorbar()
    # mpl.tight_layout()
    mpl.title('eta at t = ' + str(t/My) + ' My')
    ax.set_aspect('equal', adjustable='box')
    mpl.show()

    file.close()  
