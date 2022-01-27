# MDOODZ 7.0 makefile

# Compiler
CC = gcc

#---------------------------------------------------#
# C flags
CFLAGS = -std=c99 -D _UMFPACK_ -I ./

ifeq ($(OPT),yes)
	CFLAGS += -O3 -ftree-vectorize -funroll-loops
else
	CFLAGS += -g -Wall  -O0 -fno-inline -fno-omit-frame-pointer
endif

ifeq ($(OMP),yes)
	CFLAGS += -fopenmp -D _OMP_
else
	CFLAGS += -Wno-unknown-pragmas
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
all: Doodzi_$(MODEL)

FILES = Main_DOODZ.o FD_Jacobian.o RheologyParticles.o ChemicalRoutines.o ParticleReseeding.o Solvers.o Routines_BEN.o StokesRoutines.o StokesAssemblyCoupled.o StokesAssemblyDecoupled.o AdvectionRoutines.o RheologyDensity.o HDF5Output.o SparseTools.o ThermalRoutines.o ThermalSolver.o ParticleRoutines.o free_surf.o FlowLaws.o MemoryAllocFree.o InputOutput.o MiscFunctions.o GridRoutines.o


FILES += ./$(MODEL_PATH)/set_$(MODEL).o

Doodzi_$(MODEL): ${FILES}
	$(CC) ${FILES} -o Doodzi_$(MODEL) ${LIBS}

clean:
	rm -rf *o ./$(MODEL_PATH)/*.o Doodzi_$(MODEL)

test:
	cd ../ && python3 -m unittest test.cases.ShearTemplate
