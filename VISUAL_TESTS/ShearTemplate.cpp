#include "HDF5pp.h"
#include <Eigen/Core>
#include <iostream>
#include <cstdio>
#include "visual-tests.h"

namespace fs = std::filesystem;

using namespace Eigen;


void PlotShearTemplateReference() {
  H5p::File file = H5p::File("ShearTemplateReference.h5", "r");
  std::vector<double> params = file.read<std::vector<double>>("/Model/Params");
  const int nx = (int) params[3];
  const int nz = (int) params[4];

  std::vector<float> pressure = file.read<std::vector<float>>("/Centers/P");
  Map<MatrixXf> pressureMatrix(pressure.data(), nx - 1, nz - 1);

  std::vector<float> xcCoord = file.read<std::vector<float>>("/Model/xc_coord");
  std::vector<float> zcCoord = file.read<std::vector<float>>("/Model/zc_coord");
  std::vector<float> vxNodes = file.read<std::vector<float>>("/VxNodes/Vx");
  Map<MatrixXf> vxNodesMatrix(vxNodes.data(), nx, nz + 1);
  MatrixXf vxNodesMatrixAvg = (vxNodesMatrix.block<100, 100>(0, 1) + vxNodesMatrix.block<100, 100>(1, 1)) * 0.5;
  std::vector<float> vzNodes = file.read<std::vector<float>>("/VzNodes/Vz");
  Map<MatrixXf> vzNodesMatrix(vzNodes.data(), nx + 1, nz);
  MatrixXf vzNodesMatrixAvg = (vzNodesMatrix.block<100, 100>(1, 0) + vzNodesMatrix.block<100, 100>(1, 1)) * 0.5;

  std::vector<float> Txx = file.read<std::vector<float>>("/Centers/sxxd");
  Map<MatrixXf> TxxMatrix(Txx.data(), nx - 1, nz - 1);
  std::vector<float> Tzz = file.read<std::vector<float>>("/Centers/szzd");
  Map<MatrixXf> TzzMatrix(Tzz.data(), nx - 1, nz - 1);
  std::vector<float> Txz = file.read<std::vector<float>>("/Vertices/sxz");
  Map<MatrixXf> TxzMatrix(Txz.data(), nx, nz);

  MatrixXf Txzc = (TxzMatrix.block<100, 100>(0, 0) + TxzMatrix.block<100, 100>(1, 1) + TxzMatrix.block<100, 100>(0, 1) + TxzMatrix.block<100, 100>(1, 0)) * 0.25;
  MatrixXf Tiic = (0.5 * (TzzMatrix.array().square() + TzzMatrix.array().square()) + Txzc.array().square()).sqrt();

  std::ofstream myfile;
  myfile.open("ShearTemplateReference.dat");
  int stepsToSkip = 9;
  int j_count = 0;
  int i_count = 0;
  for (int j = 0; j < nz - 1; j++) {
    j_count++;
    for (int i = 0; i < nx - 1; i++) {
      i_count++;
      myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << pressureMatrix(i, j);
      myfile << '\t' << Tiic(i, j);
      if (j_count == 1 && i_count == 1) {
        myfile << '\t' << vxNodesMatrixAvg(i, j) << '\t' << vzNodesMatrixAvg(i, j) << std::endl;
      } else {
        myfile << std::endl;
      }
      if (i_count == stepsToSkip) {
        i_count = 0;
      }
    }
    if (j_count == stepsToSkip) {
      j_count = 0;
    }
  }
  myfile.close();

  std::FILE *GNUplotPipe = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/ShearTemplateReference.png'; data='ShearTemplateReference.dat'\" ShearTemplate.gnu", "w");
  std::fflush(GNUplotPipe);
}

void PlotShearTemplate() {
  H5p::File file = H5p::File("ShearTemplate.gzip.h5", "r");
  std::vector<double> params = file.read<std::vector<double>>("/Model/Params");
  const int nx = (int) params[3];
  const int nz = (int) params[4];

  std::vector<float> pressure = file.read<std::vector<float>>("/Centers/P");
  Map<MatrixXf> pressureMatrix(pressure.data(), nx - 1, nz - 1);

  std::vector<float> xcCoord = file.read<std::vector<float>>("/Model/xc_coord");
  std::vector<float> zcCoord = file.read<std::vector<float>>("/Model/zc_coord");
  std::vector<float> vxNodes = file.read<std::vector<float>>("/VxNodes/Vx");
  Map<MatrixXf> vxNodesMatrix(vxNodes.data(), nx, nz + 1);
  MatrixXf vxNodesMatrixAvg = (vxNodesMatrix.block<100, 100>(0, 1) + vxNodesMatrix.block<100, 100>(1, 1)) * 0.5;
  std::vector<float> vzNodes = file.read<std::vector<float>>("/VzNodes/Vz");
  Map<MatrixXf> vzNodesMatrix(vzNodes.data(), nx + 1, nz);
  MatrixXf vzNodesMatrixAvg = (vzNodesMatrix.block<100, 100>(1, 0) + vzNodesMatrix.block<100, 100>(1, 1)) * 0.5;

  std::vector<float> Txx = file.read<std::vector<float>>("/Centers/sxxd");
  Map<MatrixXf> TxxMatrix(Txx.data(), nx - 1, nz - 1);
  std::vector<float> Tzz = file.read<std::vector<float>>("/Centers/szzd");
  Map<MatrixXf> TzzMatrix(Tzz.data(), nx - 1, nz - 1);
  std::vector<float> Txz = file.read<std::vector<float>>("/Vertices/sxz");
  Map<MatrixXf> TxzMatrix(Txz.data(), nx, nz);

  MatrixXf Txzc = (TxzMatrix.block<100, 100>(0, 0) + TxzMatrix.block<100, 100>(1, 1) + TxzMatrix.block<100, 100>(0, 1) + TxzMatrix.block<100, 100>(1, 0)) * 0.25;
  MatrixXf Tiic = (0.5 * (TzzMatrix.array().square() + TzzMatrix.array().square()) + Txzc.array().square()).sqrt();

  std::ofstream myfile;
  myfile.open("ShearTemplate.dat");
  int stepsToSkip = 9;
  int j_count = 0;
  int i_count = 0;
  for (int j = 0; j < nz - 1; j++) {
    j_count++;
    for (int i = 0; i < nx - 1; i++) {
      i_count++;
      myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << pressureMatrix(i, j);
      myfile << '\t' << Tiic(i, j);
      if (j_count == 1 && i_count == 1) {
        myfile << '\t' << vxNodesMatrixAvg(i, j) << '\t' << vzNodesMatrixAvg(i, j) << std::endl;
      } else {
        myfile << std::endl;
      }
      if (i_count == stepsToSkip) {
        i_count = 0;
      }
    }
    if (j_count == stepsToSkip) {
      j_count = 0;
    }
  }
  myfile.close();

  std::FILE *GNUplotPipe = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/ShearTemplate.png'; data='ShearTemplate.dat';\" ShearTemplate.gnu", "w");
  std::fflush(GNUplotPipe);
}


void PlotShearTemplate1Reference() {
  H5p::File file = H5p::File("ShearTemplate1Reference.h5", "r");
  std::vector<double> params = file.read<std::vector<double>>("/Model/Params");
  const int nx = (int) params[3];
  const int nz = (int) params[4];

  std::vector<float> pressure = file.read<std::vector<float>>("/Centers/P");
  Map<MatrixXf> pressureMatrix(pressure.data(), nx - 1, nz - 1);

  std::vector<float> xcCoord = file.read<std::vector<float>>("/Model/xc_coord");
  std::vector<float> zcCoord = file.read<std::vector<float>>("/Model/zc_coord");
  std::vector<float> vxNodes = file.read<std::vector<float>>("/VxNodes/Vx");
  Map<MatrixXf> vxNodesMatrix(vxNodes.data(), nx, nz + 1);
  MatrixXf vxNodesMatrixAvg = (vxNodesMatrix.block<100, 100>(0, 1) + vxNodesMatrix.block<100, 100>(1, 1)) * 0.5;
  std::vector<float> vzNodes = file.read<std::vector<float>>("/VzNodes/Vz");
  Map<MatrixXf> vzNodesMatrix(vzNodes.data(), nx + 1, nz);
  MatrixXf vzNodesMatrixAvg = (vzNodesMatrix.block<100, 100>(1, 0) + vzNodesMatrix.block<100, 100>(1, 1)) * 0.5;

  std::vector<float> Txx = file.read<std::vector<float>>("/Centers/sxxd");
  Map<MatrixXf> TxxMatrix(Txx.data(), nx - 1, nz - 1);
  std::vector<float> Tzz = file.read<std::vector<float>>("/Centers/szzd");
  Map<MatrixXf> TzzMatrix(Tzz.data(), nx - 1, nz - 1);
  std::vector<float> Txz = file.read<std::vector<float>>("/Vertices/sxz");
  Map<MatrixXf> TxzMatrix(Txz.data(), nx, nz);

  MatrixXf Txzc = (TxzMatrix.block<100, 100>(0, 0) + TxzMatrix.block<100, 100>(1, 1) + TxzMatrix.block<100, 100>(0, 1) + TxzMatrix.block<100, 100>(1, 0)) * 0.25;
  MatrixXf Tiic = (0.5 * (TzzMatrix.array().square() + TzzMatrix.array().square()) + Txzc.array().square()).sqrt();

  std::ofstream myfile;
  myfile.open("ShearTemplate1Reference.dat");
  int stepsToSkip = 9;
  int j_count = 0;
  int i_count = 0;
  for (int j = 0; j < nz - 1; j++) {
    j_count++;
    for (int i = 0; i < nx - 1; i++) {
      i_count++;
      myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << pressureMatrix(i, j);
      myfile << '\t' << Tiic(i, j);
      if (j_count == 1 && i_count == 1) {
        myfile << '\t' << vxNodesMatrixAvg(i, j) << '\t' << vzNodesMatrixAvg(i, j) << std::endl;
      } else {
        myfile << std::endl;
      }
      if (i_count == stepsToSkip) {
        i_count = 0;
      }
    }
    if (j_count == stepsToSkip) {
      j_count = 0;
    }
  }
  myfile.close();

  std::FILE *GNUplotPipe = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/ShearTemplate1Reference.png'; data='ShearTemplate1Reference.dat'\" ShearTemplate.gnu", "w");
  std::fflush(GNUplotPipe);
}

void PlotShearTemplate1() {
  H5p::File file = H5p::File("ShearTemplate1.gzip.h5", "r");
  std::vector<double> params = file.read<std::vector<double>>("/Model/Params");
  const int nx = (int) params[3];
  const int nz = (int) params[4];

  std::vector<float> pressure = file.read<std::vector<float>>("/Centers/P");
  Map<MatrixXf> pressureMatrix(pressure.data(), nx - 1, nz - 1);

  std::vector<float> xcCoord = file.read<std::vector<float>>("/Model/xc_coord");
  std::vector<float> zcCoord = file.read<std::vector<float>>("/Model/zc_coord");
  std::vector<float> vxNodes = file.read<std::vector<float>>("/VxNodes/Vx");
  Map<MatrixXf> vxNodesMatrix(vxNodes.data(), nx, nz + 1);
  MatrixXf vxNodesMatrixAvg = (vxNodesMatrix.block<100, 100>(0, 1) + vxNodesMatrix.block<100, 100>(1, 1)) * 0.5;
  std::vector<float> vzNodes = file.read<std::vector<float>>("/VzNodes/Vz");
  Map<MatrixXf> vzNodesMatrix(vzNodes.data(), nx + 1, nz);
  MatrixXf vzNodesMatrixAvg = (vzNodesMatrix.block<100, 100>(1, 0) + vzNodesMatrix.block<100, 100>(1, 1)) * 0.5;

  std::vector<float> Txx = file.read<std::vector<float>>("/Centers/sxxd");
  Map<MatrixXf> TxxMatrix(Txx.data(), nx - 1, nz - 1);
  std::vector<float> Tzz = file.read<std::vector<float>>("/Centers/szzd");
  Map<MatrixXf> TzzMatrix(Tzz.data(), nx - 1, nz - 1);
  std::vector<float> Txz = file.read<std::vector<float>>("/Vertices/sxz");
  Map<MatrixXf> TxzMatrix(Txz.data(), nx, nz);

  MatrixXf Txzc = (TxzMatrix.block<100, 100>(0, 0) + TxzMatrix.block<100, 100>(1, 1) + TxzMatrix.block<100, 100>(0, 1) + TxzMatrix.block<100, 100>(1, 0)) * 0.25;
  MatrixXf Tiic = (0.5 * (TzzMatrix.array().square() + TzzMatrix.array().square()) + Txzc.array().square()).sqrt();

  std::ofstream myfile;
  myfile.open("ShearTemplate1.dat");
  int stepsToSkip = 9;
  int j_count = 0;
  int i_count = 0;
  for (int j = 0; j < nz - 1; j++) {
    j_count++;
    for (int i = 0; i < nx - 1; i++) {
      i_count++;
      myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << pressureMatrix(i, j);
      myfile << '\t' << Tiic(i, j);
      if (j_count == 1 && i_count == 1) {
        myfile << '\t' << vxNodesMatrixAvg(i, j) << '\t' << vzNodesMatrixAvg(i, j) << std::endl;
      } else {
        myfile << std::endl;
      }
      if (i_count == stepsToSkip) {
        i_count = 0;
      }
    }
    if (j_count == stepsToSkip) {
      j_count = 0;
    }
  }
  myfile.close();

  std::FILE *GNUplotPipe = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/ShearTemplate1.png'; data='ShearTemplate1.dat'\" ShearTemplate.gnu", "w");
  std::fflush(GNUplotPipe);
}