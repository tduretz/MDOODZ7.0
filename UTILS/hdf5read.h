#ifndef MDOODZ_HDF5READ_H
#define MDOODZ_HDF5READ_H
#include "hdf5.h"

typedef struct Hdf5File {
    char *fileName;
    hid_t *file_p;
    int nx;
    int nz;
} Hdf5File;

Hdf5File ReadHdf5(char *filename);

float *GetPressureArray(Hdf5File hdf5File);

#endif
