#include "HDF5pp.h"
#include <Eigen/Core>
#include <iostream>
#include <cstdio>
#include "visual-tests.h"
#include <set>
#include <filesystem>

using namespace Eigen;
using namespace std;
namespace fs = filesystem;

void BuildChart(set<fs::path> outputFiles, char* fileSuffix) {
  ofstream myfile;
  const double secondsInYear = 31556952.0;
  char datFileName[20];
  snprintf(datFileName, sizeof(datFileName), "gse%s.dat", fileSuffix);
  myfile.open(datFileName);
  double initialSize;
  for (fs::path filePath : outputFiles) {
    myfile << '\n' << '\n';
    H5p::File file = H5p::File(filePath, "r");
    vector<double> params = file.read<vector<double>>("/Model/Params");
    const double time = round(((double) params[0] / (secondsInYear * 1000))) / 1000;
    const int nx = (int) params[3];
    const int nz = (int) params[4];

    vector<float> grainSize = file.read<vector<float>>("/Centers/d");
    Map<MatrixXf> grainSizeMatrix(grainSize.data(), nx - 1, nz - 1);
    vector<float> xcCoord = file.read<vector<float>>("/Model/xc_coord");
    vector<float> zcCoord = file.read<vector<float>>("/Model/zc_coord");

    if (!initialSize) {
      initialSize = xcCoord[0];
    }
    const double extensionPercent = round((xcCoord[0] / initialSize - 1) * 100);
    myfile << "\"{/:Italic t} = " << time << " [Myr], Extension: " << extensionPercent << " % \"" << endl;
    for (int j = 0; j < nz - 1; j++) {
      for (int i = 0; i < nx - 1; i++) {
        myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << log10(grainSizeMatrix(i, j) * 1e6) << endl;
      }
    }
  }

  myfile.close();
  char command[500];
  snprintf(command, sizeof(command), R"(gnuplot -e "filename='../VISUAL_TESTS/img/gse%s.gif'" -e "data='gse%s.dat'" gse.gnu)", fileSuffix, fileSuffix);
  cout << command;
  FILE *GNUplotPipe = popen(command, "w");
  fflush(GNUplotPipe);
}

void PlotGSE() {
  std::string path = "./GSE";
  set<fs::path> outputFiles;
  for (const auto & entry : fs::directory_iterator(path)) {
    std::string filePath = entry.path();
    if (filePath.find("GSE00") == std::string::npos) {
      continue;
    }
    string path_string{entry.path().string()};
    if (path_string.find(".gzip.h5") == string::npos) {
      continue;
    }
    outputFiles.insert(entry.path());
  }
  BuildChart(outputFiles, "");
}

void PlotGSERef() {
  std::string path = "./GSE_REF";
  set<fs::path> outputFiles;
  for (const auto & entry : fs::directory_iterator(path)) {
    outputFiles.insert(entry.path());
  }
  BuildChart(outputFiles, "ref");
}