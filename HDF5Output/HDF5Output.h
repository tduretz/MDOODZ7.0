#ifndef MDOODZ_MODEL_H
#define MDOODZ_MODEL_H

typedef struct Params {
  double time;
  double Lx;
  double Lz;
  double Nx;
  double Nz;
  double dx;
  double dz;
  double dt;
} Params;

typedef struct Model_t {
  Params Params;
} Model;

typedef struct Vertices_t {
  float rho_s[ELEMENTS];
  float eta_s;
  float sxz;
  float exz;
  double rho;
  double eta;
} Vertices;

typedef struct Centers_t {

} Centers;

typedef struct VxNodes_t {

} VxNodes;

typedef struct VzNodes_t {

} VzNodes;

typedef struct Particles_t {

} Particles;

typedef struct VizGrid_t {

} VizGrid;

typedef struct Topo_t {

} Topo;

typedef struct Topo_ini_t {

} Topo_ini;

typedef struct Output_t {
  Model Model;
  Vertices Vertices;
  Centers Centers;
  VxNodes VxNodes;
  VzNodes VzNodes;
  Particles Particles;
  VizGrid VizGrid;
  Topo Topo;
  Topo_ini Topo_ini;
} Output;

#endif
