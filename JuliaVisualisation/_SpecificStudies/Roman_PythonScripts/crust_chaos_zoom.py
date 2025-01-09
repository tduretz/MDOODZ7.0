import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
from scipy.interpolate import interp2d

# Constants
My = 1e6 * 365 * 24 * 3600
Lc = 1000.0


def extract_data(file_path, data_path):
    with h5py.File(file_path, 'r') as file:
        return file[data_path][()]

def main():
    fig, ax = plt.subplots(figsize=(12, 8))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.2)

    filename = "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_800x420\\Output07840.gzip.h5"

    file   = h5py.File(filename, 'r')
    phases = file['/VizGrid/compo_hr']
    xv_ph  = file['/VizGrid/xviz_hr']
    zv_ph  = file['/VizGrid/zviz_hr']
    T      = file['Centers/T']
    data   = file['/Model/Params']
    nx   = file['Centers/nx']
    nz   = file['Centers/nz']
    xc   = file['/Model/xc_coord']
    zc   = file['/Model/zc_coord']
    Nx  = data[3].astype(int)
    Nz  = data[4].astype(int)
    t      = data[0].astype(float)
    t_ma = round(t/My, 2)
    print("the time")
    print(t_ma)
    Ncx = Nx-1
    Ncz = Nz-1

    exxd = file['/Centers/exxd']
    ezzd = file['/Centers/ezzd']
    exz = file['/Vertices/exz']
    exxd = np.reshape(exxd,(Ncz,Ncx))
    ezzd = np.reshape(ezzd,(Ncz,Ncx))
    exz = np.reshape(exz,(Nz,Nx))
    exz_c = .25*(exz[0:-1, 0:-1 ] + exz[1:, 0:-1 ] + exz[0:-1, 1: ] + exz[1:, 1: ])
    eII = (.5*exxd**2 + .5*ezzd**2 + exz_c**2)**.5
    log_eII = np.log10(eII)

    xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
    zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
    Ncx_ph = xc_ph.shape[0]
    Ncz_ph = zc_ph.shape[0]
    phases = np.reshape(phases,(Ncz_ph, Ncx_ph))
    T   = np.reshape(T,   (Ncz,Ncx))
    T -= 273.15  # Subtract 273.15 from all temperature values
    nxc  =  np.reshape(nz, (Ncz, Ncx))
    nzc  = np.reshape(nx, (Ncz, Ncx))
    # Create 2D centroid mesh
    xc2, zc2 = np.meshgrid(xc,zc)
    xc2 /= Lc
    zc2 /= Lc

    # After reshaping phases and T
    x = np.linspace(0, 1, T.shape[1])  # Assuming linear spacing
    z = np.linspace(0, 1, T.shape[0])
    f = interp2d(x, z, T, kind='linear')

    x_new = np.linspace(0, 1, phases.shape[1])
    z_new = np.linspace(0, 1, phases.shape[0])
    T_interpolated = f(x_new, z_new)

    # Define color map for phase
    cmap = ListedColormap(['white', '#ac977c', '#e2ba92', '#78c484', '#6ba5cb'])
    bounds = [-1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Plot heatmap
    ph_masked = np.ma.masked_where(phases == -1.0, phases)
    extent = [xc_ph.min()/Lc, xc_ph.max()/Lc, zc_ph.max()/Lc, zc_ph.min()/Lc]

    # Replace -inf values in log_eII
    finite_log_eII = log_eII[np.isfinite(log_eII)]
    if len(finite_log_eII) > 0:
        min_finite_log_eII = np.min(finite_log_eII)
        log_eII[log_eII == -np.inf] = min_finite_log_eII
    else:
        # If there are no finite values, use a placeholder like -10
        log_eII[log_eII == -np.inf] = -10

    # Interpolation
    x = np.linspace(0, 1, log_eII.shape[1])
    z = np.linspace(0, 1, log_eII.shape[0])
    x_new = np.linspace(0, 1, phases.shape[1])
    z_new = np.linspace(0, 1, phases.shape[0])

    f_log_eII = interp2d(x, z, log_eII, kind='cubic')
    log_eII_interpolated = f_log_eII(x_new, z_new)

    # Check for NaNs after interpolation and replace them
    log_eII_interpolated[np.isnan(log_eII_interpolated)] = -10  # Or any neutral value


    #ax.imshow(ph_masked, cmap=cmap, norm=norm, aspect='auto', extent=extent)
    # Existing code: ax.imshow for log_eII_interpolated
    image = ax.imshow(log_eII_interpolated, cmap=plt.get_cmap('inferno'), aspect='auto',
                      alpha=1, extent=extent)

    # New code: Add a color bar
    cbar = plt.colorbar(image, ax=ax)
    cbar.set_label(r'$\log(\dot{\varepsilon}_{\mathrm{II}})$ [s$^{-1}$]', rotation=270, labelpad=15)




    T_levels = [500, 1200]

    contour = ax.contour(xc_ph/Lc, zc_ph/Lc, T_interpolated, levels=T_levels, cmap='inferno', alpha=0.5)

    # Step for visualisation of quiver arrows
    stp = 10
    stp_vertical = 10

    # Interpolate phases to match the dimensions of nxc and nzc
    x_phases = np.linspace(0, 1, phases.shape[1])
    z_phases = np.linspace(0, 1, phases.shape[0])
    f_phases = interp2d(x_phases, z_phases, phases, kind='linear')

    x_new_phases = np.linspace(0, 1, nxc.shape[1])
    z_new_phases = np.linspace(0, 1, nxc.shape[0])
    phases_interpolated = f_phases(x_new_phases, z_new_phases)

    # Mask for locations where phase is not 3 in the interpolated array
    mask = (phases_interpolated[::stp_vertical, ::stp] > -0.5) & (phases_interpolated[::stp_vertical, ::stp] < 2.5)

    # Mask the quiver arrays
    nxc_masked = np.ma.masked_where(~mask, nxc[::stp_vertical, ::stp])
    nzc_masked = np.ma.masked_where(~mask, nzc[::stp_vertical, ::stp])

    n = np.sqrt(nxc_masked ** 2 + nzc_masked ** 2)

    ax.quiver(xc2[::stp_vertical, ::stp], zc2[::stp_vertical, ::stp], nxc_masked / n, nzc_masked / n, headwidth=1.,
              headaxislength=1.,
              headlength=1.,
              scale=100,
              width=0.0005)

    # Calculate the angles in radians
    angles_rad = np.arctan2(nzc_masked / n, nxc_masked / n)

    # Convert to degrees
    angles_deg = np.degrees(angles_rad)

    ax.set_ylim(-200, 10)

    ax.set_aspect('equal')
    ax.set_xlabel('x [km]', fontsize=15)
    ax.set_ylabel('y [km]', fontsize=15)


    plt.tight_layout() # Adjust the right side to make space for row labels
    plt.show()

    # Optionally save the figure
    fig.savefig("phases.png", dpi=300)


if __name__ == "__main__":
    main()
