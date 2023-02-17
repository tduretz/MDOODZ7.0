# MDOODZ7.0

Welcome to MDOODZ7.0 public repository!
This version of MDOODZ is under construction, more testing will be progressively added...

![](/misc/images/Compression_Symmetric.gif)

## [Visual Test Results](/VISUAL_TESTS/readme.md)

# Library usage

MDOODZ Source Code stored in `MDLIB` directory and compiled as a separate library `libmdoodz` 
with the public interface `mdoodz.h`.

The public interface includes a struct that stores input parameters and setup toolchains: `MdoodzInput`.
To run the simulation `MdoodzInput` must be passed to `RunMDOODZ(MdoodzInput *input)` function

### Examples:

1) Setting models with pure or simple or shear boundary conditions: [ShearTemplate](SETS/ShearTemplate.c)
2) Viscous relaxation of a free surface: [TopoBenchCase1](SETS/TopoBenchCase1.c)
3) A continental rifting and thermal solution model: [RiftingChenin](SETS/RiftingChenin.c)
4) Grain size evolution model and necking in calcite: [PinchSwellGSE](SETS/PinchSwellGSE.c)
5) Quartz-Coesite inclusion density change in Garnet: [QuartzCoesite](SETS/QuartzCoesite.c)


## Input parameters

- `inputFileName` stores the name of the `.txt` file with input parameter key = value pairs
- `model` aggregates general input parameters from `.txt` file
- `materials` aggregates input parameters from `.txt` file concerning phase properties
- `scale` aggregates input parameters from `.txt` file concerning scaling of units
- `crazyConductivity` contains parameters for the crazy conductivity of the asthenosphere for the initialisation step

### Import files

Import files are the external files that are processed by MDOODZ and contains information
such as phase transition diagrams or particle geometry. Described by `import_files_dir` and `import_file`


### Crazy conductivity

If you wish to add crazy conductivity of the asthenosphere to the initialisation step
there is a `crazyConductivity` parameter that points to the struct that aggregates
- `phases` array of phases ids that crazy conductivity should be applied to
- `nPhases` total number of phases
- `multiplier` refers to the multiplier of the effective conductivity

If your `.txt` file shares the same name as executable, 
you could extract it with the `GetSetupFileName(nargs, args)` function on Unix systems

## Setup toolchain

Structures that aggregate pointers to the functions that will be used in a runtime as callback functions.
Some of those functions must be implemented, but others if not implemented will give a default result

### BuildInitialTopography

Aggregates pointers to functions for setting up topography chain properties. 
Must have if `model.free_surf == 1`.


- `SetSurfaceZCoord` describes an altitude in relation to the x coordinate. Default value is `1.0e3 / input->scaling.L`:  flat surface will be generated
- `SetSurfacePhase` describes a topography chain particle phase id in relation to the x coordinate. Default phase `0`

### SetParticles

Aggregates pointers to functions for setting up particle properties.
Must have.


- `SetHorizontalVelocity` describes a particle Horizontal Velocity (Vx) in relation to coordinates. Default value is `-coordinates.x * input->model.EpsBG`
- `SetVerticalVelocity` describes a particle Vertical Velocity (Vz) in relation to coordinates. Default value is `coordinates.z * input->model.EpsBG`
- `SetPhase` describes a particle phase id in relation to coordinates. Default value is `0`: model will be homogeneous
- `SetTemperature` describes a particle temperature in relation to coordinates. Default value is `273.15 / input->scaling.T`: model is 0°C
- `SetGrainSize` describes a particle grain size in relation to coordinates. Default value is `0.0`
- `SetPorosity` describes a particle grain porosity in relation to coordinates. Default value is `0.0`
- `SetDensity` describes a particle grain density in relation to coordinates. Default value is set according to the particle phase 
- `SetXComponent` describes a particle X component value in relation to coordinates. Default value is `0.0`
- `SetPressure` describes a particle pressure value in relation to coordinates. Default value is `0.0`
- `SetNoise` describes a noise property value in relation to coordinates. Default value is `0.0`

## SetBCs

Aggregates pointers to functions for setting up Boundary Conditions in a mesh grid.
Must have.

- `SetBCVx` describes the type and value of the Vx point. Must be implemented. Pre-made functions from mdoodz library can be used: `SetPureShearBCVx`, `SetSimpleShearBCVx`, `SetPureOrSimpleShearBCVx` (depends on `shear_style` input parameter)
- `SetBCVz` describes the type and value of the Vz point. Must be implemented. Pre-made functions from mdoodz library can be used: `SetPureShearBCVz`, `SetSimpleShearBCVz`, `SetPureOrSimpleShearBCVz` (depends on `shear_style` input parameter)
- `SetBCPType` describes the type of the Pressure Boundary conditions point. Default one is `-1`
- `SetBCT` describes the Temperature Boundary type and value. Must be implemented if `model.isthermal == 1`
- `SetBCTNew` describes the Temperature Boundary type and value on 1d boundary array. Must be implemented if `model.isthermal == 1`. Will be deprecated

# How to build and run MDOODZ7.0?

 There are two ways to build and run MDOODZ7.0:
 - The first possibility is to use the CMake project tools. This is recommended for most users. This allows to: (1) automatically download necessary libraries (HDF5 and SuiteSparse), (2) build MDOODZ7.0 and (3) run MDOODZ7.0. 
 - The second possibility is to use the standard make project as for previous versions of MDOODZ. It is still available and recommended for users who want to manually install necesary libraries HDF5 and SuiteSparse) and control the type/versions of compilers.  

# Building and running MDOODZ7.0 using CMake

Project is ready to be built in CMake

## Prequisites and setup

In order to build MDOODZ with CMake you have to install `cmake 3.16` or newer version.
`blas`, `zlib` and `lapack`  libraries are CMake compatible and does not require any additional setup other than installing them with your package manager:

If you want to use your fixed environmental variables, set them up in a `env.cmake` file. Just copy `env.cmake.example` 
without `example` suffix.

```bash
sudo apt-get install libblas-dev liblapack-dev zlib1g-dev libhdf5-serial-dev
```

There are two ways how to use SuiteSparse:

1) Build SuiteSparse to the project:

```bash
make install-suitesparse
```

2) Specify path to the build SuiteSparse on your machine in the `env.cmake` file:

```code
set(SuiteSparse_DIR /usr/lib/x86_64-linux-gnu/)
```

If there is a SuiteSparse in a `deps` folder, it will be automatically linked with the MDOODZ. 
Otherwise it will try to find library and headers in a specified path.

## How to use

You can specify C Compiler in the `env.cmake` file. By default, it's gcc:

```code
set(C_COMPILER gcc)
```

### Add your set to the CMake

1) In a `SETS` directory you need to create executable `YourSetName.c` file and the `YourSetName.txt` setup file. Both files should have a same name.
2) In `SETS/CMakeLists.txt` file add a line with the command `add_set(YourSetName)` and it will be built as executable by default. 
3) Alternatively you could build your set with specifying it as a command line argument. Examples are given below 


### Build and run

CMake gives us a framework for building MDOODZ in developer mode or as a separate package and CTest can be used as test suite.

Makefile related to cmake is located in a root directory

To build library, tests and executables with OpenMP and in optimised mode:

```bash
make build SET=RiftingChenin
```
An alternative `.txt` input file located in the `SETS/` folder can be used if specified as `TXT`. If not stated, the default `.txt` input file will be taken.

```bash
make build SET=RiftingChenin TXT=RiftingChenin_alternative.txt
```

To explicitly set OPT (optimisation) and OMP (OpenMP). If not stated, it's OFF by default

```bash
make build-dev OMP=ON OPT=ON SET=RiftingChenin
```

Build files will be located at the `cmake-build/` directory.
Executables will be stored in a `cmake-exec/YourSetName/` directory

After building, you could run MDOODZ with:
```bash
make run SET=YourSetName
```

or run autotests:
```bash
make run-tests
```

Whole process in cluster mode with optimisation and OpenMP:

```bash
make clean build run SET=YourSetName
```

All output files generated by the runtime (`.gzip.h5`, `.dat`) will be stored in a same directory as the executable.

Setup `.txt` file will be automatically detected if it has a same name as the executable. But alternatively you can specify it explicitly:

```bash
cd cmake-exec/YourSetName
./YourSetName /path/to/setup.txt
```

# Windows 

## Setting up vcpkg

1) Install vcpkg https://vcpkg.io/en/index.html
2) By running vcpkg you will need to install HDF5 and SuiteSparse libraries by typing in 
`vcpkg install hdf5:x64-windows-static` and `vcpkg install suitesparse:x64-windows-x64-static` in a terminal.
3) Copy env.cmake.example to the same folder but without .example. Make sure that you have a correct path to vcpkg.cmake in a `CMAKE_TOOLCHAIN_FILE` param

# Building and running MDOODZ7.0 using make

It is possible to build MDOODZ7.0 using the [`makefile`](https://github.com/tduretz/MDOODZ7.0/blob/debug-strain-rate-pipo/MDLIB/makefile) which is located in the `MDLIB` folder. The `make` build does not install the necessary libraries, it is thus mandatory to manually install [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) and [HDF5](https://www.hdfgroup.org/solutions/hdf5/) prior to building MDOODZ7.0 with make. These libraries are also available via package managers (e.g., Homebrew or MacPorts). For a successful build procedure, links to compiler, library and header files should be added to the current the environement (e.g. add to `PATH`, `LIBRARY_PATH`, `C_INCLUDE_PATH`...). The excutable file will be located directly in the `MDLIB` folder and the related `YourSetName.txt` fille will be automatically copied there.

This type of build is (so far) only supported on Linux and Mac OS. 

```bash
make clean all SET=YourSetName
```
Will clean and make a fresh build without any specific optmisation.

```bash
make all SET=YourSetName
```
Will update the build (compiles only modified files) without any specific optmisation.

```bash
make clean all SET=YourSetName OPT=yes OMP=yes
```
Will clean and make a fresh build with optmisation level 3 and including OpenMP parallelism.

```bash
make all SET=YourSetName OPT=yes OMP=yes
```
Will update the build (compiles only modified files) with optmisation level 3 and including OpenMP parallelism.

