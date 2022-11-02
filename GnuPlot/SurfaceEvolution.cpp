#include "HDF5pp.h"
#include <iostream>
#include <cstdio>
#include "math.h"
#include <filesystem>

namespace fs = std::filesystem;

static std::string path = "/home/roman/CLionProjects/MDOODZ7.0/cmake-exec/RiftingChenin/topoCheck";

void PlotSurfaceEvolution() {
  std::ofstream myfile;
  myfile.open("SurfaceEvolution.dat");
  const double secondsInYear = 31556952.0;
  for (const auto & entry : fs::directory_iterator(path)) {
    std::string filePath = entry.path();
    H5p::File file = H5p::File(filePath, "r");
    std::vector<double> params = file.read<std::vector<double>>("/Model/Params");
    std::vector<float> z_mark = file.read<std::vector<float>>("/Topo/z_mark");
    float surfaceSum = 0.0;
    for (int i = 0 ; i < z_mark.size() ; i++) {
      surfaceSum += z_mark[i];
    }
    const double time = round(((double) params[0] / (secondsInYear * 1000)));
    myfile << time << '\t' << round(surfaceSum / z_mark.size()) << std::endl;
  }
  myfile.close();

  std::FILE *GNUplotPipe = popen("gnuplot SurfaceEvolution.gnu", "w");
  std::fflush(GNUplotPipe);
}

int main() {
  PlotSurfaceEvolution();
}