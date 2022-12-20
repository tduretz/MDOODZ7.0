# MDOODZ 6.0 makefile
SHELL := /bin/bash

# Compiler
#CC = /usr/local/Cellar/gcc/11.2.0_3/bin/gcc-11
CC =/opt/local/bin/gcc-mp-8

# Path to setup file folder
MODEL_PATH = 

#---------------------------------------------------#
# C flags
CFLAGS = -std=c99 -D _UMFPACK_ -I ./include/ -I ./

ifeq ($(OPT),yes)	
	CFLAGS += -O3 -ftree-vectorize -funroll-loops -finline -fomit-frame-pointer -march=native -fopenmp -D _OMP_
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
LIBS = -lz -lhdf5 -lgomp

# Related to SuiteSparse
LIBS += -L /usr/lib/ -llapack -lcxsparse -lumfpack -lcblas -lamd -lcholmod -lcolamd -lbtf -lsuitesparseconfig -lblas

# Link to openmp
ifeq ($(OMP),yes)
    LIBS += -lgomp
endif

#---------------------------------------------------#
# Rules
FILES = test_anisotropy_v1.o 


all: Doodzi_$(MODEL) 
Doodzi_$(MODEL): ${FILES}
	$(CC) ${FILES} -o test_anisotropy ${LIBS}

clean:
	rm -rf *o $(MODEL) 
