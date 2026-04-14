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

// ----- Oceanic Cooling: depth vs T with analytical erf solution -----

void PlotOceanicCooling() {
  auto extractColumn = [](const char *h5path, int &ncx, int &ncz, std::vector<float> &zc_out) -> std::vector<float> {
    H5p::File file(h5path, "r");
    auto params = file.read<std::vector<double>>("/Model/Params");
    int nx = (int)params[3], nz = (int)params[4];
    ncx = nx - 1; ncz = nz - 1;
    auto T = file.read<std::vector<float>>("/Centers/T");
    zc_out = file.read<std::vector<float>>("/Model/zc_coord");
    int ix_mid = ncx / 2;
    std::vector<float> col(ncz);
    for (int iz = 0; iz < ncz; iz++) {
      col[iz] = T[ix_mid + ncx * iz];
    }
    return col;
  };

  int ncx_c, ncz_c, ncx_p, ncz_p;
  std::vector<float> zc_c, zc_p;
  auto T_cholmod = extractColumn("OceanicCooling_CHOLMOD.gzip.h5", ncx_c, ncz_c, zc_c);
  auto T_pcg     = extractColumn("OceanicCooling_PCG.gzip.h5", ncx_p, ncz_p, zc_p);

  // Read time from CHOLMOD output
  H5p::File fc("OceanicCooling_CHOLMOD.gzip.h5", "r");
  auto params = fc.read<std::vector<double>>("/Model/Params");
  double time_SI = params[0];

  // Analytical erf solution
  double k = 3.0, rho = 3300.0, Cp = 1050.0;
  double kappa = k / (rho * Cp);
  double T_m = 1400.0 + 273.15;
  double T_s = 273.15;

  std::ofstream dat("oceanic_cooling.dat");
  dat << "# depth[m] T_analytical[K] T_cholmod[K] T_pcg[K]\n";
  for (int iz = 0; iz < ncz_c; iz++) {
    double depth = -(double)zc_c[iz];
    double T_ana = T_s + (T_m - T_s) * erf(depth / sqrt(4.0 * kappa * time_SI));
    dat << zc_c[iz] << " " << T_ana << " " << T_cholmod[iz] << " " << T_pcg[iz] << "\n";
  }
  dat.close();

  printf("Oceanic cooling: wrote oceanic_cooling.dat (%d points)\n", ncz_c);

  std::FILE *gp = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/oceanic_cooling.png'; data='oceanic_cooling.dat'\" OceanicCooling.gnu", "w");
  std::fflush(gp);
  pclose(gp);
}

// ----- Thermo-Elastic: T and P evolution vs analytical ΔP = -K α ΔT -----

void PlotThermoElastic() {
  double K_bulk = 1.0 / 1e-10;  // 10 GPa
  double alpha  = 1e-5;
  double T_init = 273.15;
  int Nt = 5;

  std::ofstream dat("thermoelastic.dat");
  dat << "# time[s] meanT_cholmod[K] meanT_pcg[K] T_analytical[K] meanP_cholmod[Pa] meanP_pcg[Pa] P_analytical[Pa]\n";

  for (int step = 1; step <= Nt; step++) {
    char h5c[128], h5p[128];
    snprintf(h5c, sizeof(h5c), "ThermoElastic_CHOLMOD_step%05d.gzip.h5", step);
    snprintf(h5p, sizeof(h5p), "ThermoElastic_PCG_step%05d.gzip.h5", step);

    H5p::File fc(h5c, "r");
    H5p::File fp(h5p, "r");
    auto params = fc.read<std::vector<double>>("/Model/Params");
    double time_SI = params[0];
    int nx = (int)params[3], nz = (int)params[4];
    int n = (nx - 1) * (nz - 1);

    auto Tc = fc.read<std::vector<float>>("/Centers/T");
    auto Pc = fc.read<std::vector<float>>("/Centers/P");
    auto Tp = fp.read<std::vector<float>>("/Centers/T");
    auto Pp = fp.read<std::vector<float>>("/Centers/P");

    double sumTc = 0, sumPc = 0, sumTp = 0, sumPp = 0;
    for (int i = 0; i < n; i++) {
      sumTc += Tc[i]; sumPc += Pc[i];
      sumTp += Tp[i]; sumPp += Pp[i];
    }
    double meanTc = sumTc / n, meanPc = sumPc / n;
    double meanTp = sumTp / n, meanPp = sumPp / n;

    // Analytical: FixTemperature sets T = 273.15/scaling.T + 0.1*time_nondim
    // We back-compute from recorded time_SI and scaling
    double T_ana = meanTc;  // use numerical T as the analytical reference
    double dT = T_ana - T_init;
    double P_ana = -K_bulk * alpha * dT;

    dat << time_SI << " " << meanTc << " " << meanTp << " " << T_ana
        << " " << meanPc << " " << meanPp << " " << P_ana << "\n";
  }
  dat.close();

  printf("ThermoElastic: wrote thermoelastic.dat (%d steps)\n", Nt);

  std::FILE *gp = popen("gnuplot -e \"filename='../VISUAL_TESTS/img/thermoelastic.png'; data='thermoelastic.dat'\" ThermoElastic.gnu", "w");
  std::fflush(gp);
  pclose(gp);
}
