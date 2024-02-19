import h5py
import numpy as np
import matplotlib.pyplot as plt

My = 1e6 * 3600 * 24 * 365.25
last_step = 2958
base_path = "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview_part\\"
particle_file_pattern = "Particles{:05d}.gzip.h5"
x_min, x_max, z_min, z_max = -101000, -100000, -18000, -17000


def read_particle_data(filename, indices=None):
    with h5py.File(filename, 'r') as file:
        dataset_length = file['Particles/x'].shape[0]

        # Filter indices to ensure they are within bounds, if indices are provided
        if indices is not None:
            valid_indices = indices[indices < dataset_length]
            x = np.array(file['Particles/x'][valid_indices])
            z = np.array(file['Particles/z'][valid_indices])
            nx = np.array(file['Particles/nx'][valid_indices])
            nz = np.array(file['Particles/nz'][valid_indices])
        else:
            x = np.array(file['Particles/x'][:])
            z = np.array(file['Particles/z'][:])
            nx = np.array(file['Particles/nx'][:])
            nz = np.array(file['Particles/nz'][:])

        t = file['Model/Params'][0].astype(float)
        t_ma = str(round(t/My, 2))
        
    return x, z, nx, nz, t_ma


def calculate_angles(nx, nz):
    n = np.sqrt(nz ** 2 + nz ** 2)

    # Calculate the angles in radians
    angles_rad = np.arctan2(nz / n, nx / n)

    # Convert to degrees
    angles_deg = np.degrees(angles_rad)
    angles_deg[angles_deg < 90] += 180
    return angles_deg


def get_tracked_particles_indices():
    filename = base_path + particle_file_pattern.format(last_step)
    x, z, nx, nz, t_ma = read_particle_data(filename, None)
    initial_indices = (x > x_min) & (x < x_max) & (z > z_min) & (z < z_max)
    indices = np.where(initial_indices)[0]
    return indices


def main():
    indices = get_tracked_particles_indices()
    steps = range(0, last_step, 2)

    positions = []
    mean_positions = []
    angles = []
    all_angles = []

    for step in steps:
        print("reading step: " + str(step))
        filename = base_path + particle_file_pattern.format(step)
        x, z, nx, nz, t_ma = read_particle_data(filename, indices)

        angles_deg = calculate_angles(nx, nz)

        mean_positions.append((np.mean(x), np.mean(z)))
        positions.append((x, z))
        angles.append(angles_deg)
        all_angles.append(angles_deg)

    # Assuming all_angles is now a list where each item is an array of angles for each step
    # Plotting
    plt.subplots_adjust(bottom=0.80, left=0.50, top=0.99)
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
    axs[1].set_ylim([0, 360])

    # For each particle, plot its angle change across steps
    for particle_index in range(len(all_angles[0])):
        particle_angles = [step_angles[particle_index] for step_angles in all_angles]
        axs[1].plot(steps, particle_angles, '-o', markersize=4)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
