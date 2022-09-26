#include "HDF5pp.h"
#include <Eigen/Core>
#include <iostream>
#include <cstdio>
#include "visual-tests.h"
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

using namespace Eigen;

void PlotRiftingCheninReference() {
  H5p::File           file    = H5p::File("RiftingCheninReference.h5", "r");
  std::vector<double> params  = file.read<std::vector<double>>("/Model/Params");
  const int           nx      = (int) params[3];
  const int           nz      = (int) params[4];
  std::vector<float>  xcCoord = file.read<std::vector<float>>("/Model/xc_coord");
  std::vector<float>  zcCoord = file.read<std::vector<float>>("/Model/zc_coord");

  std::vector<float>  exz     = file.read<std::vector<float>>("/Vertices/exz");
  Map<MatrixXf>       exzMatrix(exz.data(), nx, nz);
  std::vector<float>  exxd = file.read<std::vector<float>>("/Centers/exxd");
  Map<MatrixXf>       exxdMatrix(exxd.data(), nx - 1, nz - 1);
  MatrixXf            eII  = (0.5 * (2 * exxdMatrix.array().square() + 0.5 * (exzMatrix.block<149, 99>(1, 1).array().square() + exzMatrix.block<149, 99>(0, 0).array().square() + exzMatrix.block<149, 99>(1, 0).array().square() + exzMatrix.block<149, 99>(0, 1).array().square()))).sqrt();

  std::vector<float>  sxxd = file.read<std::vector<float>>("/Centers/sxxd");
  Map<MatrixXf>       sxxdMatrix(sxxd.data(), nx - 1, nz - 1);
  std::vector<float>  sxz = file.read<std::vector<float>>("/Vertices/sxz");
  Map<MatrixXf>       sxzMatrix(sxz.data(), nx, nz);
  MatrixXf            sII = (0.5 * (2 * sxxdMatrix.array().square() + 0.5 * (sxzMatrix.block<149, 99>(1, 1).array().square() + sxzMatrix.block<149, 99>(0, 0).array().square() + sxzMatrix.block<149, 99>(1, 0).array().square() + sxzMatrix.block<149, 99>(0, 1).array().square()))).sqrt();

  std::ofstream myfile;
  myfile.open("RiftingCheninRef.dat");
  for (int j = 0; j < nz - 1; j++) {
    for (int i = 0; i < nx - 1; i++) {
      myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << log10(eII(i, j)) << '\t' << log10(sII(i, j)) << std::endl;
    }
  }
  myfile.close();

  std::FILE *GNUplotPipe = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/RiftingCheninReference.png'\" -e \"data='RiftingCheninRef.dat'\" RiftingChenin.gnu", "w");
  std::fflush(GNUplotPipe);
}


void PlotRiftingChenin() {
  H5p::File           file    = H5p::File("RiftingChenin50.gzip.h5", "r");
  std::vector<double> params  = file.read<std::vector<double>>("/Model/Params");
  const int           nx      = (int) params[3];
  const int           nz      = (int) params[4];
  std::vector<float>  xcCoord = file.read<std::vector<float>>("/Model/xc_coord");
  std::vector<float>  zcCoord = file.read<std::vector<float>>("/Model/zc_coord");

  std::vector<float>  exz     = file.read<std::vector<float>>("/Vertices/exz");
  Map<MatrixXf>       exzMatrix(exz.data(), nx, nz);
  std::vector<float>  exxd = file.read<std::vector<float>>("/Centers/exxd");
  Map<MatrixXf>       exxdMatrix(exxd.data(), nx - 1, nz - 1);
  MatrixXf            eII  = (0.5 * (2 * exxdMatrix.array().square() + 0.5 * (exzMatrix.block<149, 99>(1, 1).array().square() + exzMatrix.block<149, 99>(0, 0).array().square() + exzMatrix.block<149, 99>(1, 0).array().square() + exzMatrix.block<149, 99>(0, 1).array().square()))).sqrt();

  std::vector<float>  sxxd = file.read<std::vector<float>>("/Centers/sxxd");
  Map<MatrixXf>       sxxdMatrix(sxxd.data(), nx - 1, nz - 1);
  std::vector<float>  sxz = file.read<std::vector<float>>("/Vertices/sxz");
  Map<MatrixXf>       sxzMatrix(sxz.data(), nx, nz);
  MatrixXf            sII = (0.5 * (2 * sxxdMatrix.array().square() + 0.5 * (sxzMatrix.block<149, 99>(1, 1).array().square() + sxzMatrix.block<149, 99>(0, 0).array().square() + sxzMatrix.block<149, 99>(1, 0).array().square() + sxzMatrix.block<149, 99>(0, 1).array().square()))).sqrt();

  std::ofstream myfile;
  myfile.open("RiftingChenin.dat");
  for (int j = 0; j < nz - 1; j++) {
    for (int i = 0; i < nx - 1; i++) {
      myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << log10(eII(i, j)) << '\t' << log10(sII(i, j)) << std::endl;
    }
  }
  myfile.close();

  std::FILE *GNUplotPipe = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/RiftingChenin.png'\" -e \"data='RiftingChenin.dat'\" RiftingChenin.gnu", "w");
  std::fflush(GNUplotPipe);
}