import h5py
import numpy as np
import matplotlib.pyplot as plt

def main():
    plt.subplots_adjust(bottom=0.80, left=0.50, top=0.99)

    filename = "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview_part\\Particles00010.gzip.h5"

    with h5py.File(filename, 'r') as file:
        # Extract the datasets and convert them to numpy arrays
        x = np.array(file['Particles/x'][:])
        z = np.array(file['Particles/z'][:])

        nx = np.array(file['Particles/nx'][:])
        nz = np.array(file['Particles/nz'][:])

        non_zero_count_nx = np.count_nonzero(nx)

        # Now you can safely perform mathematical operations
        n = np.sqrt(nx ** 2 + nz ** 2)
        nx = nx / n
        nz = nz / n

        print(1)

        # Continue with your processing...

if __name__ == "__main__":
    main()
