import h5py
import os
import numpy as np
import matplotlib.pyplot as plt

# Constants
My = 1e6 * 365 * 24 * 3600
Lc = 1000.0

output_file_pattern = "Output{:05d}.gzip.h5"
first_step = 0
step_range = 10
last_step = 7000
steps = range(first_step, last_step, step_range)
time_limit = 20  # Ma

full_lithosphere = True


def extract_data(file_path, data_path):
    with h5py.File(file_path, 'r') as file:
        return file[data_path][()]

data_paths = [
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso1\\", "label": "Anisotropy factor 1, fabric angle 80"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_0\\", "label": "Anisotropy factor 3, fabric angle 90"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_90\\", "label": "Anisotropy factor, fabric angle 0"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_800x420\\", "label": "Anisotropy factor 3 on crust and mantle, fabric angle 80"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_crust\\", "label": "Anisotropy factor 3 on crust, fabric angle 80"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_mantle\\", "label": "Anisotropy factor 3 on mantle, fabric angle 80"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_20\\", "label": "Anisotropy factor 3, fabric angle 70"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_50\\", "label": "Anisotropy factor 3, fabric angle 40"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso3_chaos1\\", "label": "Anisotropy factor 3, fabric angle random"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso1\\", "label": "Anisotropy factor 1, fabric angle 80"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso2\\", "label": "Anisotropy factor 2, fabric angle 80"},
    {"base_path": "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview3\\aniso6\\", "label": "Anisotropy factor 6, fabric angle 80"},
]

def extract_data(file_path, data_path):
    with h5py.File(file_path, 'r') as file:
        return file[data_path][()]


def generate_plot(ax, base_path, label):
    time_mas = []
    zeros_per_column_kms = []
    xc_phs = []

    for step in steps:
        filename = os.path.join(base_path, output_file_pattern.format(step))
        with h5py.File(filename, 'r') as file:
            phases = file['/VizGrid/compo_hr']
            xv_ph = file['/VizGrid/xviz_hr']
            zv_ph = file['/VizGrid/zviz_hr']
            data = file['/Model/Params']
            t = data[0].astype(float)
            t_ma = round(t / My, 2)

            if t_ma > time_limit:
                print(filename)
                break

            xc_ph = 0.5 * (xv_ph[:-1] + xv_ph[1:])
            zc_ph = 0.5 * (zv_ph[:-1] + zv_ph[1:])
            phases = np.reshape(phases, (zc_ph.shape[0], xc_ph.shape[0]))

            step_size_km = np.abs(zv_ph[1] - zv_ph[0]) / 1000.0
            if full_lithosphere:
                zeros_per_column = np.count_nonzero(np.isin(phases, [0, 1, 2]), axis=0)
            else:
                zeros_per_column = np.count_nonzero(phases == 0, axis=0)

            zeros_per_column_km = zeros_per_column * step_size_km

            time_mas.append(t_ma)
            zeros_per_column_kms.append(zeros_per_column_km)
            xc_phs.append(xc_ph / 1000)

    data_grid = np.array(zeros_per_column_kms)

    # Plot setup
    cmap = plt.get_cmap('viridis')
    mesh = ax.pcolormesh(xc_phs[0], time_mas, data_grid, cmap=cmap, shading='auto', vmin=5, vmax=140)
    ax.set_title(label)
    ax.set_xlabel('x [km]')
    ax.set_ylabel('t [Ma]')
    ax.set_ylim(0, time_limit)

    if full_lithosphere:
        plt.colorbar(mesh, ax=ax, label='Lithosphere thickness [km]')
    else:
        plt.colorbar(mesh, ax=ax, label='Crustal thickness [km]')


def main():
    num_plots = len(data_paths)
    num_cols = min(3, num_plots)
    num_rows = (num_plots + num_cols - 1) // num_cols  # Calculate the number of rows needed

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows), squeeze=False)

    for idx, data_path in enumerate(data_paths):
        ax = axes[idx // num_cols, idx % num_cols]
        generate_plot(ax, data_path['base_path'], data_path['label'])

    plt.tight_layout()
    plt.savefig('1.png')
    plt.show()


if __name__ == "__main__":
    main()
