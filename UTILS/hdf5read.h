#ifndef MDOODZ_HDF5READ_H
#define MDOODZ_HDF5READ_H
#include "hdf5.h"

typedef struct Hdf5File {
    char *fileName;
    hid_t *filePtr;
    int nx;
    int nz;
} Hdf5File;

Hdf5File ReadHdf5(char *filename);

typedef struct Coordinates {
  float *xcCoordArray;
  float *zcCoordArray;
} Coordinates;

typedef struct Nodes {
  float *VxArray;
  float *VzArray;
} Nodes;

Hdf5File ReadHdf5(char *filename);

float *GetPressureArray(Hdf5File hdf5File);

Coordinates GetCoordinates(Hdf5File hdf5File);

Nodes GetNodes(Hdf5File hdf5File);

#endif
