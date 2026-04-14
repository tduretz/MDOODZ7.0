#include "HDF5pp.h"
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include "visual-tests.h"

using namespace Eigen;

// ----- Thermal Accuracy: bar chart of L2 errors (CHOLMOD vs PCG, Gaussian + Geotherm) -----

static double computeRelL2(const std::vector<float> &num, const std::vector<double> &ana) {
  double sumDiff2 = 0.0, sumAna2 = 0.0;
  for (size_t i = 0; i < num.size(); i++) {
    double d = (double)num[i] - ana[i];
    sumDiff2 += d * d;
    sumAna2  += ana[i] * ana[i];
  }
  return (sumAna2 > 1e-30) ? sqrt(sumDiff2 / sumAna2) : sqrt(sumDiff2 / (double)num.size());
}

void PlotThermalAccuracy() {
  // --- Gaussian diffusion analytic solution ---
  auto readGaussianL2 = [](const char *h5path) -> double {
    H5p::File file(h5path, "r");
    auto params = file.read<std::vector<double>>("/Model/Params");
    double time_SI = params[0];
    int nx = (int)params[3];
    int nz = (int)params[4];
    int ncx = nx - 1, ncz = nz - 1;

    auto T = file.read<std::vector<float>>("/Centers/T");
    auto xc = file.read<std::vector<float>>("/Model/xc_coord");
    auto zc = file.read<std::vector<float>>("/Model/zc_coord");

    double k = 3.0, rho = 3300.0, Cp = 1050.0;
    double kappa = k / (rho * Cp);
    double T_bg = 500.0 + 273.15, dT = 200.0, R = 10000.0;
    double R2_4kt = R * R + 4.0 * kappa * time_SI;

    std::vector<double> T_ana(T.size());
    for (int iz = 0; iz < ncz; iz++) {
      for (int ix = 0; ix < ncx; ix++) {
        double r2 = (double)xc[ix] * xc[ix] + (double)zc[iz] * zc[iz];
        T_ana[ix + ncx * iz] = T_bg + dT * (R * R / R2_4kt) * exp(-r2 / R2_4kt);
      }
    }
    return computeRelL2(T, T_ana);
  };

  // --- Geotherm analytic solution ---
  auto readGeothermL2 = [](const char *h5path) -> double {
    H5p::File file(h5path, "r");
    auto params = file.read<std::vector<double>>("/Model/Params");
    int nx = (int)params[3];
    int nz = (int)params[4];
    double Lz = params[2];
    int ncx = nx - 1, ncz = nz - 1;

    auto T = file.read<std::vector<float>>("/Centers/T");
    auto zc = file.read<std::vector<float>>("/Model/zc_coord");

    double z_top = Lz / 2.0;
    double T_top = 273.15, T_bot = 1600.0;
    std::vector<double> T_ana(T.size());
    for (int iz = 0; iz < ncz; iz++) {
      double Tz = T_top + (T_bot - T_top) * (z_top - (double)zc[iz]) / Lz;
      for (int ix = 0; ix < ncx; ix++) {
        T_ana[ix + ncx * iz] = Tz;
      }
    }
    return computeRelL2(T, T_ana);
  };

  double L2_gauss_cholmod = readGaussianL2("GaussianDiffusionL2_CHOLMOD.gzip.h5");
  double L2_gauss_pcg     = readGaussianL2("GaussianDiffusionL2_PCG.gzip.h5");
  double L2_geotherm_cholmod = readGeothermL2("SteadyStateGeotherm_CHOLMOD.gzip.h5");
  double L2_geotherm_pcg     = readGeothermL2("SteadyStateGeotherm_PCG.gzip.h5");

  printf("Thermal accuracy L2 errors:\n");
  printf("  Gaussian  CHOLMOD: %e  PCG: %e\n", L2_gauss_cholmod, L2_gauss_pcg);
  printf("  Geotherm  CHOLMOD: %e  PCG: %e\n", L2_geotherm_cholmod, L2_geotherm_pcg);

  std::ofstream dat("thermal_accuracy.dat");
  dat << "# Test Solver L2_error\n";
  dat << "\"Gaussian\" \"CHOLMOD\" " << L2_gauss_cholmod << "\n";
  dat << "\"Gaussian\" \"PCG\" " << L2_gauss_pcg << "\n";
  dat << "\"Geotherm\" \"CHOLMOD\" " << L2_geotherm_cholmod << "\n";
  dat << "\"Geotherm\" \"PCG\" " << L2_geotherm_pcg << "\n";
  dat.close();

  std::FILE *gp = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/thermal_accuracy.png'; data='thermal_accuracy.dat'\" ThermalAccuracy.gnu", "w");
  std::fflush(gp);
  pclose(gp);
}

// ----- PCG Convergence: residual curves for cold start vs warm start -----

void PlotPCGConvergence() {
  // Read CSV files: pcg_residuals_step00000.csv (cold start) and pcg_residuals_step00004.csv (warm start)
  auto readCSV = [](const char *path, std::vector<int> &iters, std::vector<double> &abs_res, std::vector<double> &rel_res) {
    std::ifstream f(path);
    std::string line;
    std::getline(f, line);  // skip header
    while (std::getline(f, line)) {
      std::istringstream ss(line);
      std::string tok;
      std::getline(ss, tok, ','); int it = std::stoi(tok);
      std::getline(ss, tok, ','); double ar = std::stod(tok);
      std::getline(ss, tok, ','); double rr = std::stod(tok);
      iters.push_back(it);
      abs_res.push_back(ar);
      rel_res.push_back(rr);
    }
  };

  std::vector<int> it0, it4;
  std::vector<double> abs0, rel0, abs4, rel4;
  readCSV("PCGConvergence/pcg_residuals_step00000.csv", it0, abs0, rel0);
  readCSV("PCGConvergence/pcg_residuals_step00004.csv", it4, abs4, rel4);

  printf("PCG convergence: step 0 = %zu iters, step 4 = %zu iters\n", it0.size() - 1, it4.size() - 1);

  std::ofstream dat("pcg_convergence.dat");
  dat << "# iteration step0_rel step4_rel\n";
  size_t maxIts = std::max(it0.size(), it4.size());
  for (size_t i = 0; i < maxIts; i++) {
    dat << i;
    dat << " " << (i < rel0.size() ? rel0[i] : rel0.back());
    dat << " " << (i < rel4.size() ? rel4[i] : rel4.back());
    dat << "\n";
  }
  dat.close();

  std::FILE *gp = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/pcg_convergence.png'; data='pcg_convergence.dat'\" PCGConvergence.gnu", "w");
  std::fflush(gp);
  pclose(gp);
}

// ----- Interp Equivalence: 2D heatmap of |T_mode0 - T_mode3| -----

void PlotInterpEquivalence() {
  H5p::File f0("InterpEquiv_Mode0.gzip.h5", "r");
  H5p::File f3("InterpEquiv_Mode3.gzip.h5", "r");

  auto params = f0.read<std::vector<double>>("/Model/Params");
  int nx = (int)params[3], nz = (int)params[4];
  int ncx = nx - 1, ncz = nz - 1;

  auto T0 = f0.read<std::vector<float>>("/Centers/T");
  auto T3 = f3.read<std::vector<float>>("/Centers/T");
  auto xc = f0.read<std::vector<float>>("/Model/xc_coord");
  auto zc = f0.read<std::vector<float>>("/Model/zc_coord");

  double maxDiff = 0.0;
  std::ofstream dat("interp_diff.dat");
  dat << "# x z |T_mode0 - T_mode3|\n";
  for (int iz = 0; iz < ncz; iz++) {
    for (int ix = 0; ix < ncx; ix++) {
      double diff = fabs((double)T0[ix + ncx * iz] - (double)T3[ix + ncx * iz]);
      if (diff > maxDiff) maxDiff = diff;
      dat << xc[ix] << " " << zc[iz] << " " << diff << "\n";
    }
    dat << "\n";  // gnuplot pm3d block separator
  }
  dat.close();

  printf("Interp equivalence: max|T_mode0 - T_mode3| = %e\n", maxDiff);

  std::FILE *gp = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/interp_diff.png'; data='interp_diff.dat'\" InterpEquivalence.gnu", "w");
  std::fflush(gp);
  pclose(gp);
}
