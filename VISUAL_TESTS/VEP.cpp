#include "HDF5pp.h"
#include <Eigen/Core>
#include <iostream>
#include <cstdio>
#include "visual-tests.h"
#include <set>

using namespace Eigen;
using namespace std;
namespace fs = filesystem;

void BuildVEPChart(set<fs::path> outputFiles, char* fileSuffix) {
  ofstream myfile;
  const double secondsInYear = 31556952.0;
  char datFileName[20];
  snprintf(datFileName, sizeof(datFileName), "VEP%s.dat", fileSuffix);
  myfile.open(datFileName);
  for (fs::path filePath : outputFiles) {
    myfile << '\n' << '\n';
    H5p::File file = H5p::File(filePath, "r");
    vector<double> params = file.read<vector<double>>("/Model/Params");
    const double time = round((double) params[0] / secondsInYear);
    const int nx = (int) params[3];
    const int nz = (int) params[4];
    vector<float> grainSize = file.read<vector<float>>("/Centers/P");
    Map<MatrixXf> grainSizeMatrix(grainSize.data(), nx - 1, nz - 1);
    vector<float> xcCoord = file.read<vector<float>>("/Model/xc_coord");
    vector<float> zcCoord = file.read<vector<float>>("/Model/zc_coord");
    myfile << "\"{/:Italic t} = " << time << " [kyr], Extension: " << endl;
    for (int j = 0; j < nz - 1; j++) {
      for (int i = 0; i < nx - 1; i++) {
        myfile << xcCoord[i] << '\t' << zcCoord[j] << '\t' << grainSizeMatrix(i, j) / 1e6 << endl;
      }
    }
  }

  myfile.close();
  char command[500];
  snprintf(command, sizeof(command), R"(gnuplot -e "filename='../VISUAL_TESTS/img/VEP%s.gif'" -e "data='VEP%s.dat'" VEP.gnu)", fileSuffix, fileSuffix);
  cout << command;
  FILE *GNUplotPipe = popen(command, "w");
  fflush(GNUplotPipe);
}

void PlotVEP() {
  std::string path = "./VEP";
  set<fs::path> outputFiles;
  for (const auto & entry : fs::directory_iterator(path)) {
    std::string filePath = entry.path();
    if (filePath.find("VEP00") == std::string::npos) {
      continue;
    }
    string path_string{entry.path().string()};
    if (path_string.find(".gzip.h5") == string::npos) {
      continue;
    }
    outputFiles.insert(entry.path());
  }
  BuildVEPChart(outputFiles, "");
}

void PlotVEPRef() {
  std::string path = "./VEP_REF";
  set<fs::path> outputFiles;
  for (const auto & entry : fs::directory_iterator(path)) {
    outputFiles.insert(entry.path());
  }
  BuildVEPChart(outputFiles, "ref");
}