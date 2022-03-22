#ifndef MDOODZ_MATRICES_H
#define MDOODZ_MATRICES_H

typedef struct Reshape {
  int rows;
  int cols;
  int xMinLimit;
  int xMaxLimit;
  int yMinLimit;
  int yMaxLimit;
} Reshape;

Reshape InitReshape(int rows, int cols);

int *ReshapeIntArray(const int *array, Reshape *reshape);

float *ReshapeFloatArray(const float *array, Reshape *reshape);

double *ReshapeDoubleArray(const double *array, Reshape *reshape);

typedef struct Sum {
  int rows;
  int cols;
  double multiplier;
} Sum;

Sum InitSum(int rows, int cols);

int *SumIntMatrices(const int *matrix, const int *matrix2, Sum *sum);

double *SumDoubleMatrices(const double *matrix, const double *matrix2, Sum *sum);

#endif
