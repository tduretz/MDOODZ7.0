#include "HDF5pp.h"
#include <Eigen/Core>
#include <iostream>
#include <cstdio>
#include "visual-tests.h"
#include <filesystem>


using namespace Eigen;
using namespace std;

namespace fs = filesystem;

void PlotRiftingCheninReference() {
  H5p::File           file    = H5p::File("RiftingCheninReference.h5", "r");
  vector<double> params  = file.read<vector<double>>("/Model/Params");
  const int           nx      = (int) params[3];
  const int           nz      = (int) params[4];
  vector<float>  xcCoord = file.read<vector<float>>("/Model/xc_coord");
  vector<float>  zcCoord = file.read<vector<float>>("/Model/zc_coord");

  vector<float>  exz     = file.read<vector<float>>("/Vertices/exz");
  Map<MatrixXf>       exzMatrix(exz.data(), nx, nz);
  vector<float>  exxd = file.read<vector<float>>("/Centers/exxd");
  Map<MatrixXf>       exxdMatrix(exxd.data(), nx - 1, nz - 1);
  MatrixXf            eII  = (0.5 * (2 * exxdMatrix.array().square() + 0.5 * (exzMatrix.block<149, 99>(1, 1).array().square() + exzMatrix.block<149, 99>(0, 0).array().square() + exzMatrix.block<149, 99>(1, 0).array().square() + exzMatrix.block<149, 99>(0, 1).array().square()))).sqrt();

  vector<float>  sxxd = file.read<vector<float>>("/Centers/sxxd");
  Map<MatrixXf>       sxxdMatrix(sxxd.data(), nx - 1, nz - 1);
  vector<float>  sxz = file.read<vector<float>>("/Vertices/sxz");
  Map<MatrixXf>       sxzMatrix(sxz.data(), nx, nz);
  MatrixXf            sII = (0.5 * (2 * sxxdMatrix.array().square() + 0.5 * (sxzMatrix.block<149, 99>(1, 1).array().square() + sxzMatrix.block<149, 99>(0, 0).array().square() + sxzMatrix.block<149, 99>(1, 0).array().square() + sxzMatrix.block<149, 99>(0, 1).array().square()))).sqrt();

  ofstream myfile;
  myfile.open("RiftingChenin.dat");
  for (int j = 0; j < nz - 1; j++) {
    for (int i = 0; i < nx - 1; i++) {
      myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << log10(eII(i, j)) << '\t' << log10(sII(i, j)) << endl;
    }
  }
  myfile.close();

  FILE *GNUplotPipe = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/RiftingCheninReference.png'\" -e \"data='VEP%s.dat'\" RiftingChenin.gnu", "w");
  fflush(GNUplotPipe);
}


void PlotRiftingChenin() {
  H5p::File           file    = H5p::File("RiftingChenin50.gzip.h5", "r");
  vector<double> params  = file.read<vector<double>>("/Model/Params");
  const int           nx      = (int) params[3];
  const int           nz      = (int) params[4];
  vector<float>  xcCoord = file.read<vector<float>>("/Model/xc_coord");
  vector<float>  zcCoord = file.read<vector<float>>("/Model/zc_coord");

  vector<float>  exz     = file.read<vector<float>>("/Vertices/exz");
  Map<MatrixXf>       exzMatrix(exz.data(), nx, nz);
  vector<float>  exxd = file.read<vector<float>>("/Centers/exxd");
  Map<MatrixXf>       exxdMatrix(exxd.data(), nx - 1, nz - 1);
  MatrixXf            eII  = (0.5 * (2 * exxdMatrix.array().square() + 0.5 * (exzMatrix.block<149, 99>(1, 1).array().square() + exzMatrix.block<149, 99>(0, 0).array().square() + exzMatrix.block<149, 99>(1, 0).array().square() + exzMatrix.block<149, 99>(0, 1).array().square()))).sqrt();

  vector<float>  sxxd = file.read<vector<float>>("/Centers/sxxd");
  Map<MatrixXf>       sxxdMatrix(sxxd.data(), nx - 1, nz - 1);
  vector<float>  sxz = file.read<vector<float>>("/Vertices/sxz");
  Map<MatrixXf>       sxzMatrix(sxz.data(), nx, nz);
  MatrixXf            sII = (0.5 * (2 * sxxdMatrix.array().square() + 0.5 * (sxzMatrix.block<149, 99>(1, 1).array().square() + sxzMatrix.block<149, 99>(0, 0).array().square() + sxzMatrix.block<149, 99>(1, 0).array().square() + sxzMatrix.block<149, 99>(0, 1).array().square()))).sqrt();

  ofstream myfile;
  myfile.open("RiftingChenin.dat");
  for (int j = 0; j < nz - 1; j++) {
    for (int i = 0; i < nx - 1; i++) {
      myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << log10(eII(i, j)) << '\t' << log10(sII(i, j)) << endl;
    }
  }
  myfile.close();
  
  FILE *GNUplotPipe = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/RiftingChenin.png'\" RiftingChenin.gnu", "w");
  fflush(GNUplotPipe);
}