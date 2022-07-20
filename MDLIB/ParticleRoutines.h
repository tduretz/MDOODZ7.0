#ifndef MDOODZ_PARTICLEROUTINES_H
#define MDOODZ_PARTICLEROUTINES_H

typedef enum { FROM_MAT_PROP  = 0,  // interpolate from material properties structure
               FROM_PARTICLES = 1,  // interpolate straight from the particle arrays
               ANISOTROPY1    = -1, // specific for anisotropy (see Anisotropy_v2.ipynb)
               ANISOTROPY2    = -2, // specific for anisotropy (see Anisotropy_v2.ipynb)
} INTERPOLATION_MODE;

typedef struct {
  INTERPOLATION_MODE mode;
  ETA_AVG avg;
  markers particles;
  DoodzFP *mat_prop;
  grid *mesh;
  double *NodeField;
  char *NodeType;
  int prop;
  int centroid;
} InterpolateParams;

void Interpolate(int phasesCount, markers particles, DoodzFP *mat_prop, grid *mesh, double *NodeField, char *NodeType, INTERPOLATION_MODE mode, ETA_AVG avg, int prop, int centroid);

#endif