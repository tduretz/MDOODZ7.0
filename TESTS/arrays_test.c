#include "assert.h"
#include "stdlib.h"
#include "matrices.h"


static void ShouldCutMatrices() {
  const int array[16] = {
      10, 1, 3, 3,
      4, 6, 9, 7,
      8, 9, 10, 30,
      8, 9, 10, 30};
  CutCommand cutCommand;
  cutCommand.cols = 4;
  cutCommand.rows = 4;
  cutCommand.xMaxLimit = 0;
  cutCommand.xMinLimit = 1;
  cutCommand.yMinLimit = 1;
  cutCommand.yMaxLimit = -1;
  const int(*cutMatrix)[3] = ReshapeIntArray(array, &cutCommand);
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
  CutCommand cutCommand;
  cutCommand.cols = 4;
  cutCommand.rows = 4;
  cutCommand.xMinLimit = 1;
  cutCommand.xMaxLimit = -1;
  cutCommand.yMinLimit = 1;
  cutCommand.yMaxLimit = -1;
  const int(*cutMatrix)[2] = ReshapeIntArray(array, &cutCommand);
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
  CutCommand cutCommand;
  cutCommand.cols = 4;
  cutCommand.rows = 4;
  cutCommand.xMinLimit = 2;
  cutCommand.xMaxLimit = 0;
  cutCommand.yMinLimit = 2;
  cutCommand.yMaxLimit = 0;
  const int(*cutMatrix)[2] = ReshapeIntArray(array, &cutCommand);
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
  CutCommand cutCommand;
  cutCommand.cols = 4;
  cutCommand.rows = 4;
  cutCommand.xMaxLimit = -2;
  cutCommand.yMaxLimit = -2;
  cutCommand.xMinLimit = 0;
  cutCommand.yMaxLimit = 0;
  const int(*cutMatrix)[2] = ReshapeIntArray(array, &cutCommand);
  assert(cutMatrix[0][0] == 10);
  assert(cutMatrix[0][1] == 1);
  assert(cutMatrix[1][0] == 4);
  assert(cutMatrix[1][1] == 6);
  free(cutMatrix);
}

int main() {
  ShouldCutMatrices4();
  ShouldCutMatrices();
  ShouldCutMatrices2();
  ShouldCutMatrices3();
}