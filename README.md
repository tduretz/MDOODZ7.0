# MDOODZ7.0

Welcome to MDOODZ7.0 public repository!
This version of MDOODZ is under construction, more testing will be progressively added...

[Visual Test Results](/VISUAL_TESTS/readme.md)

# Library usage

MDOODZ Source Code stored in `MDLIB` directory and compiled as a separate library `libmdoodz` 
with the public interface `mdoodz.h`.

The public interface includes a struct that stores input parameters and setup toolchains: `MdoodzInstance`.
To run the simulation `MdoodzInstance` must be passed to `RunMDOODZ(MdoodzInstance *instance)` function

### Examples:

1) Minimal mechanical model: [ShearTemplate](SETS/ShearTemplate.c)
2) A free surface model: [TopoBenchCase1](SETS/TopoBenchCase1.c)
3) A free surface and thermal solution model: [RiftingPauline](SETS/RiftingPauline.c)


## Input parameters

- `inputFileName` stores the name of the `.txt` file with input parameter key = value pairs
- `model` aggregates general input parameters from `.txt` file
- `materials` aggregates input parameters from `.txt` file concerning phase properties
- `scale` aggregates input parameters from `.txt` file concerning scaling of units

If your `.txt` file shares the same name as executable, 
you could extract it with the `GetSetupFileName(nargs, args)` function on Unix systems

## Setup toolchain

Structures that aggregate pointers to the functions that will be used in a runtime as callback functions.
Some of those functions must be implemented, but others if not implemented will give a default result

### BuildInitialTopography

Aggregates pointers to functions for setting up topography chain properties. 
Must have if `model.free_surf == 1`.


- `SetSurfaceZCoord` describes an altitude in relation to the x coordinate. Default value is `1.0e3 / instance->scaling.L`:  flat surface will be generated
- `SetSurfacePhase` describes a topography chain particle phase id in relation to the x coordinate. Default phase `0`

### SetParticles

Aggregates pointers to functions for setting up particle properties.
Must have.


- `SetHorizontalVelocity` describes a particle Horizontal Velocity (Vx) in relation to coordinates. Default value is `-coordinates.x * instance->model.EpsBG`
- `SetVerticalVelocity` describes a particle Vertical Velocity (Vz) in relation to coordinates. Default value is `coordinates.z * instance->model.EpsBG`
- `SetPhase` describes a particle phase id in relation to coordinates. Default value is `0`: model will be homogeneous
- `SetTemperature` describes a particle temperature in relation to coordinates. Default value is `273.15 / instance->scaling.T`: model is 0°C
- `SetGrainSize` describes a particle grain size in relation to coordinates. Default value is `0.0`
- `SetPorosity` describes a particle grain porosity in relation to coordinates. Default value is `0.0`
- `SetDensity` describes a particle grain density in relation to coordinates. Default value is set according to the particle phase 
- `SetXComponent` describes a particle X component value in relation to coordinates. Default value is `0.0`

## SetBCs

Aggregates pointers to functions for setting up Boundary Conditions in a mesh grid.
Must have.

- `SetBCVx` describes the type and value of the Vx point. Must be implemented
- `SetBCVz` describes the type and value of the Vz point. Must be implemented
- `SetBCPType` describes the type of the Pressure Boundary conditions point. Default one is `-1`
- `SetBCT` describes the Temperature Boundary type and value. Must be implemented if `model.isthermal == 1`
- `SetBCTNew` describes the Temperature Boundary type and value on 1d boundary array. Must be implemented if `model.isthermal == 1`. Will be deprecated

## Crazy conductivity

If you wish to add crazy conductivity of the asthenosphere to the initialisation step 
there is a `crazyConductivity` parameter that points to the struct that aggregates 
- `phases` array of phases ids that crazy conductivity should be applied to
- `nPhases` total number of phases
- `multiplier` refers to the multiplier of the effective conductivity

# CMake usage

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
make build SET=RiftingPaulnie
```

To explicitly set OPT (optimisation) and OMP (OpenMP). If not stated, it's OFF by default

```bash
make build-dev OMP=ON OPT=ON SET=RiftingPaulnie
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