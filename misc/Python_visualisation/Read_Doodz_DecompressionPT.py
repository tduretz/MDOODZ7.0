from unittest.main import MAIN_EXAMPLES
import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl

file = h5py.File('/Users/tduretz/REPO_GIT/MDOODZ7.0/SOURCE/Output00200.gzip.h5', 'r')
P_mean    = file['TimeSeries/P_mean_time']
T_mean    = file['TimeSeries/T_mean_time']
time      = file['TimeSeries/Time_time']

ky = 1e3*365.25*24*3600

# Prescribed PT path: T=(P/b)^(1/a)
a      = 6.9462     
b      = 5.5398e-20
T_aim  = (P_mean[:]/1e9/b)**(1/a) + 273.15

# Predicted pressure with time
bkg_div_rate  = 1e-14
P_ini  = 3.45e9  
beta   = 1.4286e-11 
P_pred = P_ini - bkg_div_rate*time[:]/beta

fig, axs = mpl.subplots(2)
# fig.suptitle('PressureTemplate model')
#-------------------
axs[0].plot(time[:]/ky, P_pred[:]/1e9, color='r', zorder=1)
axs[0].scatter(time[::9]/ky, P_mean[::9]/1e9, marker='+', zorder=2)
axs[0].set_title('P-t diagram')
axs[0].set_xlabel('t [ky]')
axs[0].set_ylabel('P [GPa]')
axs[0].legend(['Expected', 'MDoodz'])
#-------------------
axs[1].plot(T_aim-273.15, P_mean[:]/1e9, color='r', zorder=1)
axs[1].scatter(T_mean[::9]-273.15, P_mean[::9]/1e9, marker='+', zorder=2)
axs[1].set_title('P-T diagram')
axs[1].set_xlabel('T [C]')
axs[1].set_ylabel('P [GPa]')
axs[1].legend(['Expected', 'MDoodz'])
mpl.show()

file.close()