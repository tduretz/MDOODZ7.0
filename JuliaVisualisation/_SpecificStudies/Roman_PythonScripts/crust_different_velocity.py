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

# File paths (these need to be adjusted to your actual file locations)
file_names = [
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3\\Output01200.gzip.h5", # 10Ma
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3\\Output01850.gzip.h5", # 15Ma
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3\\Output02740.gzip.h5", # 20Ma
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_fast\\Output01100.gzip.h5",
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_fast\\Output01750.gzip.h5",
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_fast\\Output02840.gzip.h5",
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_furious\\Output00870.gzip.h5",
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_furious\\Output01400.gzip.h5",
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_furious\\Output02050.gzip.h5",
]

def extract_data(file_path, data_path):
    with h5py.File(file_path, 'r') as file:
        return file[data_path][()]


def main():
    fig, axes = plt.subplots(3, 3, figsize=(12, 8)) # Adjust the layout as needed
    plt.subplots_adjust(bottom=0.80, left=0.50, top=0.99)

    column_titles = ['Extension: 1.26 cm/year', 'Extension: 2.52 cm/year', 'Extension: 3.78 cm/year']
    row_titles = []

    for i, filename in enumerate(file_names):
        file   = h5py.File(filename, 'r')
        data   = file['/Model/Params']
        t      = data[0].astype(float)
        t_ma = round(t/My, 2)
        print(filename)
        print(t_ma)
        row_titles.append(str(t_ma))


    for i, filename in enumerate(file_names):

        file   = h5py.File(filename, 'r')
        phases = file['/VizGrid/compo_hr']
        xv_ph  = file['/VizGrid/xviz_hr']
        zv_ph  = file['/VizGrid/zviz_hr']
        T      = file['Centers/T']
        strain      = file['Centers/strain']
        data   = file['/Model/Params']
        t      = data[0].astype(float)
        t_ma = round(t/My, 2)
        Nx  = data[3].astype(int)
        Nz  = data[4].astype(int)
        nx   = file['Centers/nx']
        nz   = file['Centers/nz']
        xc   = file['/Model/xc_coord']
        zc   = file['/Model/zc_coord']
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
        strain   = np.reshape(strain,   (Ncz,Ncx))
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
        f2 = interp2d(x, z, strain, kind='linear')

        x_new = np.linspace(0, 1, phases.shape[1])
        z_new = np.linspace(0, 1, phases.shape[0])
        T_interpolated = f(x_new, z_new)
        strain_interpolated = f2(x_new, z_new)

        # Plotting
        ax = axes[i % 3, i // 3]

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


        ax.imshow(ph_masked, cmap=cmap, norm=norm, aspect='auto', extent=extent)
        ax.imshow(strain_interpolated, cmap=plt.get_cmap('inferno'), aspect='auto',
                  alpha=0, extent=extent, vmin=0, vmax=10)

        


        T_levels = [500, 1200]

        contour = ax.contour(xc_ph/Lc, zc_ph/Lc, T_interpolated, levels=T_levels, cmap='inferno', alpha=0.5)

        # Step for visualisation of quiver arrows
        stp = 9
        stp_vertical = 9

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
                  headlength=1.,)

        ax.set_ylim(-100, 3)

        row = i // 3
        col = i % 3

        ax.set_aspect('equal')

        if row == 0:
            axes[row, col].set_title(column_titles[col], fontsize=20)

            # Set x-ticks for the last row
        if col < 2:
            ax.set_xticks([])  # Remove x-ticks for all but the last row

        # Set y-ticks for the first column
        if row > 0:
            ax.set_yticks([])  # Remove y-ticks for all but the first column

        # Set x and y labels
        if i % 3 == 2:
            ax.set_xlabel('x [km]', fontsize=15)
        if i // 3 == 0:
            ax.set_ylabel('y [km]', fontsize=15)

        ax.text(1.02, 0.5, str(t_ma) + ' Ma',
                transform=ax.transAxes, ha='left', va='center',
                rotation='vertical', fontsize=15)



    plt.tight_layout(rect=[0, 0.1, 1, 1]) # Adjust the right side to make space for row labels
    plt.show()

    # Optionally save the figure
    fig.savefig("phases.png", dpi=300)


if __name__ == "__main__":
    main()
