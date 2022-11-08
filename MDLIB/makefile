# MDOODZ 6.0 makefile
SHELL := /bin/bash

# Compiler
CC = /usr/local/Cellar/gcc/11.2.0_3/bin/gcc-11

# Path to setup file folder
MODEL_PATH = ../SETS

#---------------------------------------------------#
# C flags
CFLAGS = -std=c99 -D _UMFPACK_ -I ./include/ -I ./

ifeq ($(OPT),yes)	
	CFLAGS += -O3 -ftree-vectorize -funroll-loops -finline -fomit-frame-pointer -march=native
else
	CFLAGS += -g -Wall  -O0 -fno-inline -fno-omit-frame-pointer
endif

ifeq ($(OMP),yes)
	CFLAGS += -fopenmp -D _OMP_
else
	CFLAGS += -Wno-unknown-pragmas 
endif

ifeq ($(VG),yes)
	CFLAGS += -D _VG_ -Wno-format-zero-length
endif

ifeq ($(NEW_INPUT),yes)
	CFLAGS += -D _NEW_INPUT_
endif

CFLAGS += -Wno-unused-variable  -Wno-comment
#---------------------------------------------------#
# Libraries
LIBS = -lz -lhdf5

# Related to SuiteSparse
LIBS += -L /usr/lib/ -llapack -lcxsparse -lumfpack -lcblas -lamd -lcholmod -lcolamd -lbtf -lsuitesparseconfig -lblas

# Link to openmp
ifeq ($(OMP),yes)
    LIBS += -lgomp
endif

#---------------------------------------------------#
# Rules
FILES = AnisotropyRoutines.o Main_DOODZ.o FD_Jacobian.o  RheologyParticles.o ChemicalRoutines.o ParticleReseeding.o Solvers.o Routines_BEN.o StokesRoutines.o StokesAssemblyCoupled.o StokesAssemblyDecoupled.o AdvectionRoutines.o RheologyDensity.o HDF5Output.o SparseTools.o ThermalRoutines.o ThermalSolver.o ParticleRoutines.o free_surf.o FlowLaws.o MemoryAllocFree.o InputOutput.o MiscFunctions.o GridRoutines.o setup.o $(MODEL_PATH)/$(MODEL).o 


all: Doodzi_$(MODEL) 
Doodzi_$(MODEL): ${FILES}
	$(CC) ${FILES} -o $(MODEL) ${LIBS}

clean:
	cp $(MODEL_PATH)/$(MODEL).txt . 
	rm -rf *o ./$(MODEL_PATH)/*.o $(MODEL) 
