import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl

def descend_obj(obj,sep='\t'):
    """
    Iterate through groups in a HDF5 file and prints the groups and datasets names and datasets attributes
    """
    if type(obj) in [h5py._hl.group.Group,h5py._hl.files.File]:
        for key in obj.keys():
            print(sep,'-',key,':',obj[key])
            descend_obj(obj[key],sep=sep+'\t')
    elif type(obj)==h5py._hl.dataset.Dataset:
        for key in obj.attrs.keys():
            print(sep+'\t','-',key,':',obj.attrs[key])

path="./"
ext=".gzip.h5"
for i in range(0,51,5):
    file = h5py.File(path+'Output'+"%05d" % i+ext, 'r')
    # print(file.keys())
    # print(descend_obj(file))
    P    = file['Centers/P']
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
    # P = np.log10(P)

    exxd = file['/Centers/exxd']
    ezzd = file['/Centers/ezzd']
    exz = file['/Vertices/exz']

    exxd = np.reshape(exxd,(Ncz,Ncx))
    ezzd = np.reshape(ezzd,(Ncz,Ncx))
    exz = np.reshape(exz,(Nz,Nx))
    
    exz_c = .25*(exz[0:-1, 0:-1 ] + exz[1:, 0:-1 ] + exz[0:-1, 1: ] + exz[1:, 1: ])

    eII = (.5*exxd**2 + .5*ezzd**2 + exz_c**2)**.5
    log_eII = np.log10(eII)

    mpl.clf()
    mpl.pcolor(xc,zc,log_eII, vmin=-18, vmax=-12)
    ax = mpl.gca()
    ax.set_aspect('equal')
    mpl.set_cmap('jet')
    mpl.colorbar()
    mpl.tight_layout()
    # mpl.show()
    mpl.title('eII')
    mpl.savefig(path+'eII'+str(i)+'.png')

