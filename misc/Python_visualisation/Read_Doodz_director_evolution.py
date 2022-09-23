import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
import os

dir_path = '/Users/romankulakov/CLionProjects/MDOODZ7/cmake-exec/SimpleShearAnisoHomo/'

directory = os.fsencode(dir_path)

nxs = []
nzs = []
Txys = []

files = os.listdir(directory)
files.sort()

print(files)

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if not filename.endswith(".h5"):
        continue
    file = h5py.File(f'{dir_path}{filename}', 'r')
    nx   = file['Centers/nx']
    nz   = file['Centers/nz']
    data = file['/Model/Params']
    Txx  = file['Centers/sxxd']
    Tzz  = file['Centers/szzd']
    Txz  = file['Vertices/sxz']
    Nx  = data[3].astype(int)
    Nz  = data[4].astype(int)
    Ncx = Nx-1
    Ncz = Nz-1
    Txx = np.reshape(Txx, (Ncz,Ncx)).transpose()
    Tzz = np.reshape(Tzz, (Ncz,Ncx)).transpose()
    Txz = np.reshape(Txz, (Nz ,Nx )).transpose()

    nxs.append(np.mean(nx))
    nzs.append(np.mean(nz))


    Txzc = 0.25*( Txz[0:-1,0:-1] + Txz[0:-1,1:]  + Txz[1:,0:-1]  + Txz[1:,1:] )
    Txys.append(np.mean(Txzc))


mpl.plot(nxs)
mpl.show()

