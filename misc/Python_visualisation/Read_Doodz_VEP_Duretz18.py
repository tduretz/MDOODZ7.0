import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl
from os.path import expanduser
home = expanduser("~")

# Path relative to your home directory
path = "/REPO_GIT/MDOODZ7.0_OK/SOURCE/"

# File numbers
file_start = 20
file_step  = 5
file_end   = 30

# Allocate arrays for logs
num_files  = int( (file_end-file_start)/file_step + 1)
time_all   = np.zeros((num_files), dtype=float)
max_topo   = np.zeros((num_files), dtype=float)

# Loop through files
i = 0
for step in range(file_start, file_end+1, file_step):
    filename    = home + path + 'Output%05d' % step + '.gzip.h5'
    file        = h5py.File(filename, 'r')
    rx              = file['/Iterations/rx_abs']
    LogIsNewtonStep = file['/Iterations/LogIsNewtonStep']
    NumberSteps     = file['/Iterations/NumberSteps']
    data            = file['/Model/Params']
    rx1 = rx[0:NumberSteps[0]+1]
    q   = np.log( np.abs( (rx1[3:] + rx1[2:-1]) / (rx1[2:-1]- rx1[1:-2]) )  ) / np.log( np.abs( (rx1[2:-1] + rx1[1:-2]) / (rx1[1:-2]- rx1[0:-3]) )  )
    print( q )

    a = np.log10(rx1[1:]) - np.log10(rx1[0:-1])


    # a = np.log10(rx1[0:-2]) + np.log10(rx1[2:]) - 2* np.log10(rx1[1:-1])

    print(a)

    fig, axs = mpl.subplots(1)
    axs.plot( range(0, NumberSteps[0]+1, 1), np.log10(rx[0:NumberSteps[0]+1]), color='r', zorder=1)
    # axs.set_ylabel('mean(Tii) [MPa]')
    # axs.set_xlabel('t [ky]')
    # axs.legend(['Expected', 'MDoodz'])
    axs.set_title('Step %03d'% step)

    i          += 1
    # # if step == file_end:
    # #     Time_time     = file['/TimeSeries/Time_time']
    # #     Tii_mean_time = file['/TimeSeries/Tii_mean_time']
    # #     # print(Time_time[:])
    # #     # print(Tii_mean_time[:])
    # #     # Figure
    # #     ky       = 1e3*365.25*24*3600

    # #     fig, axs = mpl.subplots(1)
    # #     axs.plot(Time_time[:]/ky, Tii_mean_time[:]/1e6, color='r', zorder=1)
    # #     axs.set_ylabel('mean(Tii) [MPa]')
    # #     axs.set_xlabel('t [ky]')
    # #     # axs.legend(['Expected', 'MDoodz'])
    # #     # axs.set_title('Case 1 test from Crameri et al. (2012)')
    # #     mpl.show()
       
    file.close()

mpl.show()


