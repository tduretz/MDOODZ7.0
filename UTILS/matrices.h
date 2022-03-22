#ifndef MDOODZ_MATRICES_H
#define MDOODZ_MATRICES_H

typedef struct CutCommand {
  int rows;
  int cols;
  int xMinLimit;
  int xMaxLimit;
  int yMinLimit;
  int yMaxLimit;
} CutCommand;

int *ReshapeIntArray(int *array, CutCommand *cutCommand);

int *ReshapeFloatArray(float *array, CutCommand *cutCommand);

int *ReshapeDoubleArray(double *array, CutCommand *cutCommand);

#endif
