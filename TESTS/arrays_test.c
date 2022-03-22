#include "assert.h"
#include "stdlib.h"
#include "matrices.h"
#include "stdio.h"


static void ShouldCutMatrices() {
  const int array[16] = {
      10, 1, 3, 3,
      4, 6, 9, 7,
      8, 9, 10, 30,
      8, 9, 10, 30};
  Reshape reshape;
  reshape.cols = 4;
  reshape.rows = 4;
  reshape.xMaxLimit = 0;
  reshape.xMinLimit = 1;
  reshape.yMinLimit = 1;
  reshape.yMaxLimit = -1;
  const int(*cutMatrix)[3] = ReshapeIntArray(array, &reshape);
  assert(cutMatrix[0][0] == 6);
  assert(cutMatrix[0][2] == 7);
  assert(cutMatrix[1][1] == 10);
  assert(cutMatrix[1][2] == 30);
  free(cutMatrix);
}

static void ShouldCutMatrices2() {
  const int array[16] = {
      10, 1, 3, 3,
      4, 6, 9, 7,
      8, 9, 10, 30,
      8, 9, 10, 30};
  Reshape reshape;
  reshape.cols = 4;
  reshape.rows = 4;
  reshape.xMinLimit = 1;
  reshape.xMaxLimit = -1;
  reshape.yMinLimit = 1;
  reshape.yMaxLimit = -1;
  const int(*cutMatrix)[2] = ReshapeIntArray(array, &reshape);
  assert(cutMatrix[0][0] == 6);
  assert(cutMatrix[0][1] == 9);
  assert(cutMatrix[1][0] == 9);
  assert(cutMatrix[1][1] == 10);
  free(cutMatrix);
}

static void ShouldCutMatrices3() {
  const int array[16] = {
      10, 1, 3, 3,
      4, 6, 9, 7,
      8, 9, 10, 30,
      8, 9, 10, 30};
  Reshape reshape;
  reshape.cols = 4;
  reshape.rows = 4;
  reshape.xMinLimit = 2;
  reshape.xMaxLimit = 0;
  reshape.yMinLimit = 2;
  reshape.yMaxLimit = 0;
  const int(*cutMatrix)[2] = ReshapeIntArray(array, &reshape);
  assert(cutMatrix[0][0] == 10);
  assert(cutMatrix[0][1] == 30);
  assert(cutMatrix[1][0] == 10);
  assert(cutMatrix[1][1] == 30);
  free(cutMatrix);
}

static void ShouldCutMatrices4() {
  const int array[16] = {
      10, 1, 3, 3,
      4, 6, 9, 7,
      8, 9, 10, 30,
      8, 9, 10, 30};
  Reshape reshape;
  reshape.cols = 4;
  reshape.rows = 4;
  reshape.xMaxLimit = -2;
  reshape.yMaxLimit = -2;
  reshape.xMinLimit = 0;
  reshape.yMaxLimit = 0;
  const int(*cutMatrix)[2] = ReshapeIntArray(array, &reshape);
  assert(cutMatrix[0][0] == 10);
  assert(cutMatrix[0][1] == 1);
  assert(cutMatrix[1][0] == 4);
  assert(cutMatrix[1][1] == 6);
  free(cutMatrix);
}

static void ShouldSumMatrices() {
  const int array[16] = {
      10, 1, 3, 3,
      4, 6, 9, 7,
      8, 9, 10, 30,
      8, 9, 10, 30};
  Reshape reshape;
  reshape.cols = 4;
  reshape.rows = 4;
  reshape.xMinLimit = 0;
  reshape.xMaxLimit = -2;
  reshape.yMaxLimit = -2;
  reshape.yMinLimit = 0;
  const int(*cutMatrix)[2] = ReshapeIntArray(array, &reshape);
  Reshape reshape2;
  reshape2.cols = 4;
  reshape2.rows = 4;
  reshape2.xMinLimit = 2;
  reshape2.xMaxLimit = 0;
  reshape2.yMinLimit = 2;
  reshape2.yMaxLimit = 0;
  const int(*cutMatrix2)[2] = ReshapeIntArray(array, &reshape2);
  Sum sum;
  sum.cols = 2;
  sum.rows = 2;
  const int(*sumMatrix)[2] = SumIntMatrices(cutMatrix, cutMatrix2, &sum);
  assert(sumMatrix[0][0] == 20);
  assert(sumMatrix[0][1] == 31);
  assert(sumMatrix[1][0] == 14);
  assert(sumMatrix[1][1] == 36);
  free(cutMatrix);
  free(cutMatrix2);
  free(sumMatrix);
}

static void ShouldSumDoubleMatrices() {
  const double array[16] = {
      10.0, 1.0, 3.0, 3.0,
      4.0, 6.0, 9.0, 7.0,
      8.0, 9.0, 10.0, 30.0,
      8.0, 9.0, 10.0, 30.0};
  Reshape reshape;
  reshape.cols = 4;
  reshape.rows = 4;
  reshape.xMinLimit = 0;
  reshape.xMaxLimit = -2;
  reshape.yMaxLimit = -2;
  reshape.yMinLimit = 0;
  const double (*cutMatrix)[2] = ReshapeDoubleArray(array, &reshape);
  printf("cutMatrix: %f %f\n", cutMatrix[0][0], cutMatrix[0][1]);
  printf("cutMatrix: %f %f\n", cutMatrix[1][0], cutMatrix[1][1]);
  Reshape reshape2;
  reshape2.cols = 4;
  reshape2.rows = 4;
  reshape2.xMinLimit = 2;
  reshape2.xMaxLimit = 0;
  reshape2.yMinLimit = 2;
  reshape2.yMaxLimit = 0;
  const double (*cutMatrix2)[2] = ReshapeDoubleArray(array, &reshape2);
  printf("cutMatrix2: %f %f\n", cutMatrix2[0][0], cutMatrix2[0][1]);
  printf("cutMatrix2: %f %f\n", cutMatrix2[1][0], cutMatrix2[1][1]);
  Sum sum;
  sum.cols = 2;
  sum.rows = 2;
  const double (*sumMatrix)[2] = SumDoubleMatrices(cutMatrix, cutMatrix2, &sum);
  printf("sumMatrix: %f %f\n", sumMatrix[0][0], sumMatrix[0][1]);
  printf("sumMatrix: %f %f\n", sumMatrix[1][0], sumMatrix[1][1]);
  assert(sumMatrix[0][0] == 20.0);
  assert(sumMatrix[0][1] == 31.0);
  assert(sumMatrix[1][0] == 14.0);
  assert(sumMatrix[1][1] == 36.0);
  free(cutMatrix);
  free(cutMatrix2);
  free(sumMatrix);
  exit(1);
}

int main() {
  ShouldSumDoubleMatrices();
}