#include "mdoodz.h"
#include "math.h"

typedef struct {
    double x, y;
} Point;

typedef struct {
    int numVertices;
    Point *vertices;
} Polygon;

// Function to check if a point is inside a polygon
bool isPointInsidePolygon(Point point, Polygon polygon) {
    int i, j;
    bool isInside = false;
    for (i = 0, j = polygon.numVertices - 1; i < polygon.numVertices; j = i++) {
        if ((polygon.vertices[i].y > point.y) != (polygon.vertices[j].y > point.y) &&
            (point.x < (polygon.vertices[j].x - polygon.vertices[i].x) * (point.y - polygon.vertices[i].y) /
                       (polygon.vertices[j].y - polygon.vertices[i].y) + polygon.vertices[i].x)) {
            isInside = !isInside;
        }
    }
    return isInside;
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {

  // Fig. 10) Garnet pure shear
  double Angle    = 60.0 / 180.0 * M_PI;

  // scaling geometry values
  // Angle doesn't need scaling (in radians btw)
  
  double A[2] = { -0.0924021300036714, -0.0665524544007828 };
  double B[2] = {  0.0198515118129514, -0.1222508262945290 };
  double C[2] = {  0.1346758477169800, -0.0125678785653086 };
  double D[2] = {  0.0172808177255486,  0.1176806218631410 };
  double E[2] = { -0.0906883339454028,  0.0508425755906477 };

  Point testPoint = {coordinates.x, coordinates.z};
  Polygon testPolygon;
  testPolygon.numVertices = 5;
  testPolygon.vertices = (Point[]){{A[0], A[1]}, {B[0], B[1]}, {C[0], C[1]}, {D[0], D[1]}, {E[0], E[1]}};

  return isPointInsidePolygon(testPoint, testPolygon);
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (zeroC) / input->scaling.T;
  if (1 == 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

double SetAnisoAngle(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double radius = input->model.user1 / input->scaling.L;

  // Fig. 10) Garnet pure shear
  double Angle    = 60.0 / 180.0 * M_PI;

  // scaling geometry values
  // Angle doesn't need scaling (in radians btw)
  
  double A[2] = { -0.0924021300036714, -0.0665524544007828 };
  double B[2] = {  0.0198515118129514, -0.1222508262945290 };
  double C[2] = {  0.1346758477169800, -0.0125678785653086 };
  double D[2] = {  0.0172808177255486,  0.1176806218631410 };
  double E[2] = { -0.0906883339454028,  0.0508425755906477 };

  Point testPoint = {coordinates.x, coordinates.z};
  Polygon testPolygon;
  testPolygon.numVertices = 5;
  testPolygon.vertices = (Point[]){{A[0], A[1]}, {B[0], B[1]}, {C[0], C[1]}, {D[0], D[1]}, {E[0], E[1]}};

  if (isPointInsidePolygon(testPoint, testPolygon)) {
    return Angle * 180.0 / M_PI + 90; // corresponding aniso_angle
  } else {
    return rand()*360; // random value in matrix
  }
}

int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
                   .SetAnisoAngle         = SetAnisoAngle,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ("Fig10_bench.txt", &setup);
}
