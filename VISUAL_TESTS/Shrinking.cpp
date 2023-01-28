#include "HDF5pp.h"
#include <Eigen/Core>
#include <iostream>
#include <cstdio>
#include "visual-tests.h"
#include <set>

namespace fs = std::filesystem;

using namespace Eigen;
using namespace std;

void BuildShrinkingChart(set<fs::path> outputFiles, char* fileSuffix) {
  ofstream myfile;
  const double secondsInYear = 31556952.0;
  char datFileName[20];
  snprintf(datFileName, sizeof(datFileName), "shrinking%s.dat", fileSuffix);
  myfile.open(datFileName);
  for (fs::path filePath : outputFiles) {
    myfile << '\n' << '\n';
    H5p::File file = H5p::File(filePath, "r");
    vector<double> params = file.read<vector<double>>("/Model/Params");
    const double time = round((params[0] / 31556952.0) * 1000) / 1000;
    const int nx = (int) params[3];
    const int nz = (int) params[4];
    
    vector<float> grainSize = file.read<vector<float>>("/Centers/X");
    Map<MatrixXf>  xMatrix(grainSize.data(), nx - 1, nz - 1);
    vector<float> xcCoord = file.read<vector<float>>("/Model/xc_coord");
    vector<float> zcCoord = file.read<vector<float>>("/Model/zc_coord");

    std::vector<float> Txx = file.read<std::vector<float>>("/Centers/sxxd");
    Map<MatrixXf> TxxMatrix(Txx.data(), nx - 1, nz - 1);
    std::vector<float> Tzz = file.read<std::vector<float>>("/Centers/szzd");
    Map<MatrixXf> TzzMatrix(Tzz.data(), nx - 1, nz - 1);
    std::vector<float> Txz = file.read<std::vector<float>>("/Vertices/sxz");
    Map<MatrixXf> TxzMatrix(Txz.data(), nx, nz);

    MatrixXf Txzc = (TxzMatrix.block<50, 50>(0, 0) + TxzMatrix.block<50, 50>(1, 1) + TxzMatrix.block<50, 50>(0, 1) + TxzMatrix.block<50, 50>(1, 0)) * 0.25;
    MatrixXf Tiic = (0.5 * (TzzMatrix.array().square() + TzzMatrix.array().square()) + Txzc.array().square()).sqrt();

    myfile << "\"{/:Italic t} = " << time << " [yr]" << endl;
    for (int j = 0; j < nz - 1; j++) {
      for (int i = 0; i < nx - 1; i++) {
        float stress = log10(Tiic(i, j));
        if (isnan(stress) || isinf(stress)) {
          stress = 0;
        }
        myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << xMatrix(i, j) << '\t' << stress << endl;
      }
    }
  }

  myfile.close();
  char command[500];
  snprintf(command, sizeof(command), R"(gnuplot -e "filename='../VISUAL_TESTS/img/shrinking%s.gif'" -e "data='shrinking%s.dat'" shrinking.gnu)", fileSuffix, fileSuffix);
  cout << command;
  FILE *GNUplotPipe = popen(command, "w");
  fflush(GNUplotPipe);
}

void PlotShrinkingGif() {
  std::string path = "./Shrinking";
  set<fs::path> outputFiles;
  for (const auto & entry : fs::directory_iterator(path)) {
    std::string filePath = entry.path();
    outputFiles.insert(entry.path());
  }
  BuildShrinkingChart(outputFiles, "");
}

void PlotShrinkingGifRef() {
  std::string path = "./Shrinking_REF";
  set<fs::path> outputFiles;
  for (const auto & entry : fs::directory_iterator(path)) {
    std::string filePath = entry.path();
    outputFiles.insert(entry.path());
  }
  BuildShrinkingChart(outputFiles, "Ref");
}