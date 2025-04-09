import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Constants
My = 1e6 * 365 * 24 * 3600
Lc = 1000.0

output_file_pattern = "Output{:05d}.gzip.h5"
first_step = 0
step_range = 10
last_step = 7000
steps = range(first_step, last_step, step_range)
time_limit = 20  # Ma


def extract_data(file_path, data_path):
    with h5py.File(file_path, 'r') as file:
        return file[data_path][()]

data_paths = [
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso6\\", "label": "Anisotropy factor 6"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_800x420\\", "label": "Anisotropy factor 3"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso1\\", "label": "Anisotropy factor 1"},
]

def extract_data(file_path, data_path):
    with h5py.File(file_path, 'r') as file:
        return file[data_path][()]


def generate_plot(ax, base_path, label):
    time_mas = []
    xc_phs = []
    velocity_means = []

    for step in steps:
        filename = os.path.join(base_path, output_file_pattern.format(step))
        with h5py.File(filename, 'r') as file:
            vx_mark = file['/Topo/Vx_mark'][:]
            x_mark = file['/Topo/x_mark'][:] / 1000  # Convert x from meters to kilometers
            data = file['/Model/Params'][:]
            t = data[0].astype(float)
            t_ma = round(t / My, 2)

            if t_ma > time_limit:
                print(f"Time limit reached at step {step}: {filename}")
                break

            time_mas.append(t_ma)
            xc_phs.append(x_mark)  # Assuming x_mark is already structured or averaged per timestep
            velocity_means.append(vx_mark * 31.536e6 * 100)  # Assuming vx_mark is averaged or processed similarly

    # Assuming time_mas, xc_phs, and velocity_means are aligned and structured correctly
    T, X = np.meshgrid(np.array(time_mas), np.array(xc_phs[0]))  # Use the first timestep as the x basis
    V = np.array(velocity_means).T  # Transpose to align dimensions correctly

    cmap = plt.get_cmap('seismic')
    mesh = ax.pcolormesh(X, T, V, cmap=cmap, shading='auto',  vmin=-1, vmax=1)
    ax.set_title(label)
    ax.set_xlabel('X Coordinate (km)')
    ax.set_ylabel('Time (Ma)')
    ax.set_ylim(0, time_limit)
    ax.set_xlim(-200, 200)
    plt.colorbar(mesh, ax=ax, label='Surface Velocity [mm/yr]')


def main():
    num_plots = len(data_paths)
    num_cols = min(3, num_plots)
    num_rows = (num_plots + num_cols - 1) // num_cols  # Calculate the number of rows needed

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows), squeeze=False)

    for idx, data_path in enumerate(data_paths):
        ax = axes[idx // num_cols, idx % num_cols]
        generate_plot(ax, data_path['base_path'], data_path['label'])

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
