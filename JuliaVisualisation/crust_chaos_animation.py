import h5py
import imageio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
from scipy.interpolate import interp2d
import os

# Constants
My = 1e6 * 365 * 24 * 3600
Lc = 1000.0

base_path = "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview_part\\aniso3_chaos1"
output_file_pattern = "Output{:05d}.gzip.h5"
particles_file_pattern = "Particles{:05d}.gzip.h5"
first_step = 0
step_range = 2
last_step = 2958
steps = range(first_step, last_step, step_range)

tracked_particle = 1193596


def extract_data(file_path, data_path):
    with h5py.File(file_path, 'r') as file:
        return file[data_path][()]


# Adjusted main function to create a GIF
def create_frames_and_gif():
    frames_dir = "/mnt/data/frames/"
    if not os.path.exists(frames_dir):
        os.makedirs(frames_dir)

    for step in steps:
        particles_path = os.path.join(base_path, particles_file_pattern.format(step))

        particles_file   = h5py.File(particles_path, 'r')
        tracked_particle_x = particles_file['Particles/x'][292039]
        tracked_particle_z = particles_file['Particles/z'][292039]

        file_path = os.path.join(base_path, output_file_pattern.format(step))
        print("reading step: " + str(step))
        if not os.path.exists(file_path):
            print(f"File {file_path} not found.")
            continue  # Skip if file doesn't exist

        file   = h5py.File(file_path, 'r')
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

        fig, ax = plt.subplots(figsize=(12, 8))


        ax.imshow(ph_masked, cmap=cmap, norm=norm, aspect='auto', extent=extent)
        ax.imshow(log_eII_interpolated, cmap=plt.get_cmap('inferno'), aspect='auto',
                  alpha=0.5, extent=extent)

        T_levels = [500, 1200]

        contour = ax.contour(xc_ph/Lc, zc_ph/Lc, T_interpolated, levels=T_levels, cmap='inferno', alpha=0.5)

        # Step for visualisation of quiver arrows
        stp = 1
        stp_vertical = 1

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



        # Calculate the angles in radians
        angles_rad = np.arctan2(nzc_masked / n, nxc_masked / n)

        # Convert to degrees
        angles_deg = np.degrees(angles_rad)

        ax.set_ylim(-100, 3)
        ax.set_xlim(-200, 200)

        ax.set_aspect('equal')
        ax.set_xlabel('x [km]', fontsize=15)
        ax.set_ylabel('y [km]', fontsize=15)

        ax.plot(tracked_particle_x / Lc, tracked_particle_z / Lc, 'ro')
        plt.title(f"Step {step} Time: {t_ma} Ma")

        frame_path = os.path.join(frames_dir, f"frame_{step:05d}.png")
        plt.savefig(frame_path)
        plt.close(fig)

    # Create GIF from frames
    frames = [imageio.imread(os.path.join(frames_dir, f)) for f in sorted(os.listdir(frames_dir)) if f.endswith('.png')]
    gif_path = "/mnt/data/animation.gif"
    imageio.mimsave(gif_path, frames, fps=10)

    return gif_path


def main():
    fig, ax = plt.subplots(figsize=(12, 8))
    plt.subplots_adjust(bottom=0.80, left=0.50, top=0.99)
    create_frames_and_gif()


if __name__ == "__main__":
    main()
