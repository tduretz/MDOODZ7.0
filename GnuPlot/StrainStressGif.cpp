#include "HDF5pp.h"
#include <Eigen/Core>
#include <iostream>
#include <cstdio>
#include <set>
#include <filesystem>

using namespace Eigen;
using namespace std;
namespace fs = filesystem;

void BuildVEPChart(set<fs::path> outputFiles) {
  ofstream myfile;
  const double secondsInYear = 31556952.0;
  myfile.open("StrainStressGif.dat");
  for (fs::path filePath : outputFiles) {
    myfile << '\n' << '\n';
    H5p::File file = H5p::File(filePath, "r");
    vector<double> params = file.read<vector<double>>("/Model/Params");
    const double time = round((double) params[0] / secondsInYear);
    int           nx      = (int) params[3];
    int           nz      = (int) params[4];
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


    myfile << "\"{/:Italic t} = " << time / 1000000 << " [Myr], Extension: " << endl;
    for (int j = 0; j < nz - 1; j++) {
      for (int i = 0; i < nx - 1; i++) {
        myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << log10(eII(i, j)) << '\t' << log10(sII(i, j)) << std::endl;
      }
    }
  }

  myfile.close();
  char command[500];
  snprintf(command, sizeof(command), R"(gnuplot -e "filename='../VISUAL_TESTS/img/StrainStressGif.gif'" -e "data='StrainStressGif.dat'" StrainStressGif.gnu)");
  cout << command;
  FILE *GNUplotPipe = popen(command, "w");
  fflush(GNUplotPipe);
}

int main() {
  std::string path = "/home/roman/CLionProjects/MDOODZ7.0/cmake-exec/CollisionIra/result6";
  set<fs::path> outputFiles;
  for (const auto & entry : fs::directory_iterator(path)) {
    outputFiles.insert(entry.path());
  }
  BuildVEPChart(outputFiles);
}