import h5py
import numpy as np
import matplotlib.pyplot as plt

My     = 1e6*3600*24*365.25


def read_particle_data(filename, indices=None):
    with h5py.File(filename, 'r') as file:
        # Extract only the data for particles of interest if indices are provided
        if indices is not None:
            x = np.array(file['Particles/x'][indices])
            z = np.array(file['Particles/z'][indices])
            nx = np.array(file['Particles/nx'][indices])
            nz = np.array(file['Particles/nz'][indices])
        else:
            x = np.array(file['Particles/x'][:])
            z = np.array(file['Particles/z'][:])
            nx = np.array(file['Particles/nx'][:])
            nz = np.array(file['Particles/nz'][:])
        t = file['Model/Params'][0].astype(int)
    return x, z, nx, nz, t


def calculate_angles(nx, nz):
    angles_rad = np.arctan2(nz, nx)
    angles_deg = np.degrees(angles_rad)
    angles_deg = np.mod(angles_deg, 360)
    angles_deg[angles_deg > 180] -= 360
    angles_deg = np.abs(angles_deg)
    return angles_deg


def main():
    # Adjust plot layout
    plt.subplots_adjust(bottom=0.80, left=0.50, top=0.99)

    # Define file pattern and range
    base_path = "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\particles\\"
    file_pattern = "Particles{:05d}.gzip.h5"
    steps = range(0, 12, 2)  # From 0 to 10 inclusive, with a step of 2

    initial_indices = None
    positions = []
    mean_positions = []
    angles = []

    for step in steps:
        filename = base_path + file_pattern.format(step)
        x, z, nx, nz, t = read_particle_data(filename, initial_indices)
        time = str(t/My)

        if initial_indices is None:
            initial_indices = (x > 5000) & (x < 6000) & (z < -3000) & (z > -4000)
            x, z, nx, nz = x[initial_indices], z[initial_indices], nx[initial_indices], nz[initial_indices]

        angles_deg = calculate_angles(nx, nz)

        mean_positions.append((np.mean(x), np.mean(z)))
        positions.append((x, z))
        angles.append(angles_deg)

    true_count = np.sum(initial_indices)

    # Plotting
    fig, axs = plt.subplots(2, 1, figsize=(10, 8))

    # Plot x vs z for selected particles
    mean_x, mean_z = zip(*mean_positions)
    axs[0].plot(mean_x, mean_z, '-o', markersize=4, label='Mean Trajectory')
    axs[0].set_title('Mean Particle Trajectory (x vs z)')
    axs[0].set_xlabel('Mean x')
    axs[0].set_ylabel('Mean z')
    axs[0].legend()

    # Plot angle changes Assuming angles is a list of numpy arrays, where each array contains angles for all tracked
    # particles at a given step
    mean_angles_per_step = [np.mean(angle_array) for angle_array in angles]

    # Correcting the plotting call
    axs[1].plot(steps, mean_angles_per_step, '-o', markersize=4)
    axs[1].set_title('Particle Angle Changes')
    axs[1].set_xlabel('Step')
    axs[1].set_ylabel('Angle (degrees)')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
