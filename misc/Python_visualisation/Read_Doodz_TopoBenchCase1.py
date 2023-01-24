import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
from os.path import expanduser
home = expanduser("~")

# Path relative to your home directory
path = "/REPO_GIT/MDOODZ7.0/MDLIB/"

# File numbers
file_start = 0
file_step  = 10
file_end   = 100

# Allocate arrays for logs
num_files  = int( (file_end-file_start)/file_step + 1)
time_all   = np.zeros((num_files), dtype=float)
max_topo   = np.zeros((num_files), dtype=float)

# Loop through files
i = 0
for step in range(file_start, file_end+1, file_step):
    filename    = home + path + 'Output%05d' % step + '.gzip.h5'
    file        = h5py.File(filename, 'r')
    z_mark      = file['/Topo/z_mark']
    data        = file['/Model/Params']
    time_all[i] = data[0]
    max_topo[i] = np.max(z_mark)
    i          += 1
    file.close()

# Analytical solution from Crameri et al. (2012)
ky       = 1e3*365.25*24*3600
hmax     = 7e3
h_dot    = 0.2139e-11
exp_topo = hmax*np.exp(-h_dot*time_all)

# Compute error
err = abs(exp_topo-max_topo)/abs(exp_topo)*100
print(err)

# Figure
fig, axs = mpl.subplots(1)
axs.plot(time_all/ky, exp_topo/1e3, color='r', zorder=1)
axs.plot(time_all/ky, max_topo/1e3, color='b', zorder=2, marker='+', linewidth=0)
axs.set_xlabel('t [ky]')
axs.set_ylabel('h [m]')
axs.legend(['Expected', 'MDoodz'])
axs.set_title('Case 1 test from Crameri et al. (2012)')
mpl.show()

