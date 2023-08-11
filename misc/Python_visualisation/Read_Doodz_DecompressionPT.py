from unittest.main import MAIN_EXAMPLES
import h5py as h5py
import numpy as np
import matplotlib.pyplot as mpl

file   = h5py.File('/Users/tduretz/REPO/MDOODZ7.0/MDLIB/Output00080.gzip.h5', 'r')
P_mean = file['TimeSeries/P_mean_time']
T_mean = file['TimeSeries/T_mean_time']
time   = file['TimeSeries/Time_time']

ky = 1e3*365.25*24*3600

# Prescribed PT path: T=(P/b)^(1/a)
a           = 1.02e-2     # fit coefficient a for T=(1/a)*ln(P/b) 
b           = 2.523e-3 # fit coefficient b for T=(1/a)*ln(P/b) 
# T_aim  = (P_mean[:]/1e9/b)**(1/a) + 273.15
T_aim  = np.log(P_mean[:]/1e9/b)*(1/a) + 273.15

# Predicted pressure with time
bkg_div_rate  = 1e-14
P_ini  = 4.3e9  
beta   = 5.8824e-12 
P_pred = P_ini - bkg_div_rate*time[:]/beta

print(P_mean[:])
print(time[:]/ky)

sp = 1

fig, axs = mpl.subplots(2)
# fig.suptitle('PressureTemplate model')
#-------------------
axs[0].plot(time[:]/ky, P_pred[:]/1e9, color='r', zorder=1)
axs[0].scatter(time[:80:sp]/ky, P_mean[:80:sp]/1e9, marker='+', zorder=2)
axs[0].set_title('P-t diagram')
axs[0].set_xlabel('t [ky]')
axs[0].set_ylabel('P [GPa]')
axs[0].legend(['Expected', 'MDoodz'])
#-------------------
axs[1].plot(T_aim-273.15, P_mean[:]/1e9, color='r', zorder=1)
axs[1].scatter(T_mean[:80:sp]-273.15, P_mean[:80:sp]/1e9, marker='+', zorder=2)
axs[1].set_title('P-T diagram')
axs[1].set_xlabel('T [C]')
axs[1].set_ylabel('P [GPa]')
axs[1].legend(['Expected', 'MDoodz'])
mpl.show()

file.close()

