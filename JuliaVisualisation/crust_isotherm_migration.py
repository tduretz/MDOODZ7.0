import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import os

# Constants
My = 1e6 * 365 * 24 * 3600
Lc = 1000.0

base_path = "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_20"
output_file_pattern = "Output{:05d}.gzip.h5"
first_step = 0
step_range = 10
last_step = 2740
steps = range(first_step, last_step, step_range)


def extract_data(file_path, data_path):
    with h5py.File(file_path, 'r') as file:
        return file[data_path][()]

def main():
    time_ma = []
    closest_x_coords = []

    for step in steps:
        filename = os.path.join(base_path, output_file_pattern.format(step))

        file   = h5py.File(filename, 'r')
        phases = file['/VizGrid/compo_hr']
        xv_ph  = file['/VizGrid/xviz_hr']
        zv_ph  = file['/VizGrid/zviz_hr']
        T      = file['Centers/T']
        data   = file['/Model/Params']
        Nx  = data[3].astype(int)
        Nz  = data[4].astype(int)
        t      = data[0].astype(float)
        t_ma = round(t/My, 2)
        Ncx = Nx-1
        Ncz = Nz-1

        xc_ph  = 0.5*(xv_ph[0:-1] + xv_ph[1:])
        zc_ph  = 0.5*(zv_ph[0:-1] + zv_ph[1:])
        Ncx_ph = xc_ph.shape[0]
        Ncz_ph = zc_ph.shape[0]
        phases = np.reshape(phases,(Ncz_ph, Ncx_ph))
        T   = np.reshape(T,   (Ncz,Ncx))
        T -= 273.15  # Subtract 273.15 from all temperature values

        # After reshaping phases and T
        x = np.linspace(0, 1, T.shape[1])  # Assuming linear spacing
        z = np.linspace(0, 1, T.shape[0])
        f = interp2d(x, z, T, kind='linear')

        x_new = np.linspace(0, 1, phases.shape[1])
        z_new = np.linspace(0, 1, phases.shape[0])
        T_interpolated = f(x_new, z_new)

        close_points_mask = (T_interpolated >= 1195) & (T_interpolated <= 1200)
        if np.any(close_points_mask):
            # Get z-indices of these points
            z_indices = np.where(close_points_mask)[0]
            # Find the maximum z-index
            max_z_index = np.max(z_indices)
            # Find all x-indices corresponding to this z-index
            x_indices_at_max_z = np.where(close_points_mask[max_z_index])[0]
            # If there are multiple, choose the first one (or any other logic you prefer)
            chosen_x_index = x_indices_at_max_z[0]
            # Convert to real x-coordinate
            max_z_x_coord = xc_ph[chosen_x_index] / Lc
            closest_x_coords.append(max_z_x_coord)
            time_ma.append(t_ma)


    # Plotting
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(closest_x_coords, time_ma)
    ax.set_xlabel('X-coordinate of max Z at T=1300 [km]', fontsize=15)
    ax.set_ylabel('Time [Ma]', fontsize=15)
    plt.tight_layout()
    plt.show()
    fig.savefig("time_vs_x_at_max_z_T_1300.png", dpi=300)


if __name__ == "__main__":
    main()
