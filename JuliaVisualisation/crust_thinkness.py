import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap
from scipy.interpolate import interp2d

# Constants
My = 1e6 * 365 * 24 * 3600
Lc = 1000.0


def extract_data(file_path, data_path):
    with h5py.File(file_path, 'r') as file:
        return file[data_path][()]

file_groups = [
    [
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso1\\Output02550.gzip.h5", "Layer angle = 10°, Isotropic: δ = 1"),
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso2\\Output02680.gzip.h5", "Layer angle = 10°, Moderate anisotropy: δ = 2"),
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso6\\Output03040.gzip.h5", "Layer angle = 10°, High anisotropy: δ = 6"),
    ],
    [
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3\\Output02710.gzip.h5", "Layer angle = 10°, Anisotropic mantle and crust: δ = 3"),
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_crust\\Output02510.gzip.h5", "Layer angle = 10°, Anisotropic crust: δ = 3"),
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_mantle\\Output02690.gzip.h5", "Layer angle = 10°, Anisotropic mantle: δ = 3"),
    ],
    [
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3\\Output02710.gzip.h5", "Layer angle = 10°: δ = 3"),
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_50\\Output02670.gzip.h5", "Layer angle = 40°: δ = 3"),
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_20\\Output02740.gzip.h5", "Layer angle = 70°: δ = 3"),
    ],
    [
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_chaos1\\Output02700.gzip.h5", "Random layer angle 1: δ = 3"),
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_chaos2\\Output02710.gzip.h5", "Random layer angle 1: δ = 3"),
        ( "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_chaos3\\Output02710.gzip.h5", "Random layer angle 1: δ = 3"),
    ],
]



def main():
    color_palettes = ['Blues', 'Oranges', 'Greens', 'Purples']
    fig, axes = plt.subplots(len(file_groups), 1, figsize=(12, 8), sharex=True, sharey=True)
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.1)

    for group_index, group in enumerate(file_groups):
        cmap = get_cmap(color_palettes[group_index])
        colors = cmap(np.linspace(0.4, 1, len(group)))  # Adjust as needed

        for file_index, (filename, legend) in enumerate(group):
            # Load the file and data as in your main function
            with h5py.File(filename, 'r') as file:
                file   = h5py.File(filename, 'r')
                phases = file['/VizGrid/compo_hr']
                xv_ph  = file['/VizGrid/xviz_hr']
                zv_ph  = file['/VizGrid/zviz_hr']
                T      = file['Centers/T']
                data   = file['/Model/Params']
                xc   = file['/Model/xc_coord']
                zc   = file['/Model/zc_coord']
                Nx  = data[3].astype(int)
                Nz  = data[4].astype(int)
                Ncx = Nx-1
                Ncz = Nz-1

                xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
                zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
                Ncx_ph = xc_ph.shape[0]
                Ncz_ph = zc_ph.shape[0]
                phases = np.reshape(phases,(Ncz_ph, Ncx_ph))
                zeros_per_column = np.count_nonzero(phases == 0, axis=0)

                # Assuming zc_ph is evenly spaced, calculate the vertical step size
                step_size_km = np.abs(zc_ph[1] - zc_ph[0]) / 1000.0  # Convert from meters to kilometers if necessary

                # Multiply zeros_per_column by the step size to convert to kilometers
                zeros_per_column_km = zeros_per_column * step_size_km

                # Assuming the rest of the processing is similar to your main function
                step_size_km = np.abs(zv_ph[1] - zv_ph[0]) / 1000.0
                zeros_per_column = np.count_nonzero(phases == 0, axis=0)
                zeros_per_column_km = zeros_per_column * step_size_km
                xc_ph = 0.5 * (xv_ph[:-1] + xv_ph[1:])

                # Plot
                ax = axes[group_index]
                ax.plot(xc_ph / 1000.0, zeros_per_column_km, color=colors[file_index], label=legend)
                ax.legend()
                ax.set_xlim(-200, 200)
                ax.set_ylim(0, 35)


    # Plotting
    ax.plot(xc_ph / 1000.0, zeros_per_column_km)  # Convert xc_ph to kilometers if necessary
    ax.set_xlabel('x [km]', fontsize=15)
    fig.text(0.01, 0.5, 'Crust thickness [km]', va='center', rotation='vertical', fontsize=15)
    plt.show()


if __name__ == "__main__":
    main()
