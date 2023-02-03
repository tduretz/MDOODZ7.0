#include "HDF5pp.h"
#include <iostream>
#include <cstdio>
#include "visual-tests.h"
#include "math.h"

namespace fs = std::filesystem;

void PlotTopoBenchCase1() {
  std::string path = ".";
  std::ofstream myfile;
  myfile.open("TopoBenchCase1.dat");
  const double secondsInYear = 31556952.0;
  for (const auto & entry : fs::directory_iterator(path)) {
    std::string filePath = entry.path();
    if (filePath.find("TopoBenchCase000") == std::string::npos) {
      continue;
    }
    H5p::File file = H5p::File(filePath, "r");
    std::vector<double> params = file.read<std::vector<double>>("/Model/Params");
    std::vector<float> z_mark = file.read<std::vector<float>>("/Topo/z_mark");
    const double time = (double) params[0] / (secondsInYear * 1000);
    myfile << time << '\t' << z_mark[0] << std::endl;
  }
  myfile.close();

  std::FILE *GNUplotPipe = popen("gnuplot TopoBenchCase1.gnu", "w");
  std::fflush(GNUplotPipe);
}