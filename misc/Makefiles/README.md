
Check the [Documentation](https://github.com/tduretz/MDOODZ6.0/blob/master/Documentation/MDOODZ_docu.pdf)

0. go to SOURCE/ folder
1. copy one makefile from Makefile/ folder into SOURCE/
2. edit the makefile for your specific computer

## Simple shear (power-law viscous)

Compile:
```bash
make clean all MODEL=Shear_pwl OPT=yes OMP=yes
```

Run:
```bash
./Doodzi_Shear_pwl Shear_pwl.txt
```

## Lithosphere deformation
Compile: make clean all MODEL=LithoScale OPT=yes OMP=yes

Run:
```bash
./Doodzi_LithoScale LithoScale.txt
```

Visualize using either Matlab, Python, Julia or whatever language that can handle HDF5 files and enjoy!

# Prequisites

So far MDOODZ has been successfully built on LINUX/UNIX and MAC OS systems. The code can be built with GCC compiler from GNU (http://gcc.gnu.org) or with ICC compiler from Intel MKL library (https://software.intel.com/en-us/intel-mkl).
The code relies on two libraries: <br>
1. SuiteSparse provides efficient linear algebra and matrix manipulation routines. It is available at: http://www.suitesparse.com <br>
2. HDF5 is the main format for output files and is readable into MATLAB. It is available at: http://www.hdfgroup.org <br>

Optionally, it is possible to use Gnuplot to monitor progress of non-linear iterations. Find it there: http://www.gnuplot.info/
