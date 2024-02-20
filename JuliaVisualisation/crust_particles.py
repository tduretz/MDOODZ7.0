import h5py
import numpy as np
import matplotlib.pyplot as plt

My = 1e6 * 3600 * 24 * 365.25
first_step = 0
step_range = 2
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
        
    return x, z, nx, nz


def read_time(filename):
    with h5py.File(filename, 'r') as file:
        t = file['Model/Params'][0].astype(float)
        t_ma = round(t/My, 2)
    return t_ma


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
    x, z, nx, nz = read_particle_data(filename, None)
    initial_indices = (x > x_min) & (x < x_max) & (z > z_min) & (z < z_max)
    indices = np.where(initial_indices)[0]
    return indices


def select_closest_to_whole_numbers(data):
    new_dict = {}
    new_dict[0] = 0.0

    current_value = 0.0

    for key, value in data.items():
        if value < current_value:
            continue
        else:
            new_dict[key] = value
            current_value += 1.0

    key, value = data.popitem()
    new_dict[key] = value

    return new_dict


def main():
    indices = get_tracked_particles_indices()
    min_value = np.min(indices)
    indices = np.array([min_value])

    steps = range(first_step, last_step, step_range)

    my_to_steps = {}

    for step in steps:
        filename = base_path + particle_file_pattern.format(step)
        t_ma = read_time(filename)
        my_to_steps[step] = t_ma

    my_to_steps = select_closest_to_whole_numbers(my_to_steps)

    positions = []
    mean_positions = []
    angles = []
    all_angles = []

    for key, value in my_to_steps.items():
        step = key
        print("reading step: " + str(step))
        filename = base_path + particle_file_pattern.format(step)
        x, z, nx, nz = read_particle_data(filename, indices)

        angles_deg = calculate_angles(nx, nz)

        mean_positions.append((np.mean(x), np.mean(z)))
        positions.append((x, z))
        angles.append(angles_deg)
        all_angles.append(angles_deg)

    # Assuming all_angles is now a list where each item is an array of angles for each step
    # Plotting
    plt.subplots_adjust(bottom=0.80, left=0.50, top=0.99)
    fig, axs = plt.subplots(2, 1, figsize=(10, 8))

    steps_for_positions = list(my_to_steps.keys())

    # Plotting mean x vs z for selected particles with annotations
    for i, (x, z) in enumerate(mean_positions):
        # Plot each point
        axs[0].plot(x, z, 'o', markersize=4)

        # Annotate each point with its corresponding time value from my_to_steps
        step = steps_for_positions[i]  # Get the step corresponding to this position
        time_ma = my_to_steps[step]  # Get the time in Ma for this step
        axs[0].annotate(f'{time_ma} Ma', (x, z), textcoords="offset points", xytext=(0,10), ha='center')

    # Set titles and labels
    axs[0].set_title('Mean Particle Trajectory (x vs z)')
    axs[0].set_xlabel('Mean x')
    axs[0].set_ylabel('Mean z')

    # Add a legend
    axs[0].legend()

    # Plot angle changes Assuming angles is a list of numpy arrays, where each array contains angles for all tracked
    # particles at a given step
    mean_angles_per_step = [np.mean(angle_array) for angle_array in angles]

    # Correcting the plotting call
    axs[1].plot(my_to_steps.keys(), mean_angles_per_step, '-o', markersize=4)
    axs[1].set_title('Particle Angle Changes')
    axs[1].set_xlabel('Step')
    axs[1].set_ylabel('Angle (degrees)')
    axs[1].set_ylim([0, 360])

    # For each particle, plot its angle change across steps
    for particle_index in range(len(all_angles[0])):
        particle_angles = [step_angles[particle_index] for step_angles in all_angles]
        axs[1].plot(my_to_steps.keys(), particle_angles, '-o', markersize=4)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
