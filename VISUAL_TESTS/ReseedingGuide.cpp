/*
 * ReseedingGuide.cpp — Visual-test runner and TSV extractor for the
 * reseeding visual guide.
 *
 * Runs rigid-body rotation advection with reseed_mode 0/1/2 and edge-case
 * configurations, then extracts particle data to TSV for gnuplot.
 */

#include "HDF5pp.h"
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <filesystem>
#include <vector>
#include <set>
#include <algorithm>
#include <sys/wait.h>
#include <unistd.h>
#include "visual-tests.h"

extern "C" {
#include "mdoodz.h"
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace fs = std::filesystem;
using namespace std;

// ============================================================
//  MDOODZ callbacks — rigid-body rotation (2 phases, disk)
// ============================================================

static const double OMEGA    = 2.0 * M_PI;
static const double DISK_CX  = 0.25;
static const double DISK_CZ  = 0.0;
static const double DISK_R   = 0.1;

static int RGSetPhase(MdoodzInput *input, Coordinates coords) {
  double x = coords.x - DISK_CX / input->scaling.L;
  double z = coords.z - DISK_CZ / input->scaling.L;
  double r = DISK_R   / input->scaling.L;
  return (x * x + z * z < r * r) ? 1 : 0;
}

static double RGSetDensity(MdoodzInput *input, Coordinates coords, int phase) {
  return input->materials.rho[phase];
}

static double RGSetHorizontalVelocity(MdoodzInput *input, Coordinates coords) {
  return -OMEGA * coords.z * input->scaling.L / input->scaling.V;
}

static double RGSetVerticalVelocity(MdoodzInput *input, Coordinates coords) {
  return OMEGA * coords.x * input->scaling.L / input->scaling.V;
}

// ============================================================
//  MDOODZ callbacks — closed-box convection roll
//  v·n = 0 on all walls, div(v) = 0
// ============================================================

static const double CONV_AMP = 200.0;

static double CBSetHorizontalVelocity(MdoodzInput *input, Coordinates coords) {
  double x = coords.x * input->scaling.L;  // dimensional
  double z = coords.z * input->scaling.L;
  double vx = CONV_AMP * sin(M_PI * (x + 0.5)) * cos(M_PI * (z + 0.5));
  return vx / input->scaling.V;
}

static double CBSetVerticalVelocity(MdoodzInput *input, Coordinates coords) {
  double x = coords.x * input->scaling.L;
  double z = coords.z * input->scaling.L;
  double vz = -CONV_AMP * cos(M_PI * (x + 0.5)) * sin(M_PI * (z + 0.5));
  return vz / input->scaling.V;
}

// --- 3-phase multiphase variant: two inclusions + matrix ---

static int RGSetPhaseMultiphase(MdoodzInput *input, Coordinates coords) {
  double x = coords.x;
  double z = coords.z;
  double r = DISK_R / input->scaling.L;
  // Disk 1 at (+0.25, 0) -> phase 1
  double x1 = x - 0.25 / input->scaling.L;
  double z1 = z - 0.0;
  if (x1 * x1 + z1 * z1 < r * r) return 1;
  // Disk 2 at (-0.25, 0) -> phase 2
  double x2 = x + 0.25 / input->scaling.L;
  double z2 = z - 0.0;
  if (x2 * x2 + z2 * z2 < r * r) return 2;
  return 0;
}

// Helper: read int8 datasets that HDF5pp doesn't natively support
static vector<int8_t> readInt8Dataset(const string &h5path, const string &dsPath) {
  H5::H5File f(h5path, H5F_ACC_RDONLY);
  H5::DataSet ds = f.openDataSet(dsPath);
  size_t n = ds.getSpace().getSelectNpoints();
  vector<int8_t> data(n);
  ds.read(data.data(), H5::PredType::NATIVE_INT8);
  return data;
}

// ============================================================
//  TSV Extractors
// ============================================================

// Extract all particles from an HDF5 file, writing TSV with:
//   x  z  phase  generation  index
static void extractParticlesForCell(const string &h5path, const string &tsvPath,
                                    double xlo, double xhi, double zlo, double zhi) {
  H5p::File file(h5path, "r");
  vector<float> px = file.read<vector<float>>("/Particles/x");
  vector<float> pz = file.read<vector<float>>("/Particles/z");
  vector<int8_t> pp = readInt8Dataset(h5path, "/Particles/phase");
  vector<int8_t> pg = readInt8Dataset(h5path, "/Particles/generation");

  ofstream out(tsvPath);
  out << "# x\tz\tphase\tgeneration\tindex\n";
  for (int i = 0; i < (int)px.size(); i++) {
    if (px[i] >= xlo && px[i] <= xhi && pz[i] >= zlo && pz[i] <= zhi) {
      out << px[i] << "\t" << pz[i] << "\t" << (int)pp[i] << "\t" << (int)pg[i] << "\t" << i << "\n";
    }
  }
  out.close();
}

// Extract all particles from an HDF5 file (no spatial filter)
static void extractAllParticles(const string &h5path, const string &tsvPath) {
  H5p::File file(h5path, "r");
  vector<float> px = file.read<vector<float>>("/Particles/x");
  vector<float> pz = file.read<vector<float>>("/Particles/z");
  vector<int8_t> pp = readInt8Dataset(h5path, "/Particles/phase");
  vector<int8_t> pg = readInt8Dataset(h5path, "/Particles/generation");

  ofstream out(tsvPath);
  out << "# x\tz\tphase\tgeneration\tindex\n";
  for (int i = 0; i < (int)px.size(); i++) {
    out << px[i] << "\t" << pz[i] << "\t" << (int)pp[i] << "\t" << (int)pg[i] << "\t" << i << "\n";
  }
  out.close();
}

// Compute per-fine-cell particle count histogram for an HDF5 snapshot.
// Bins active particles (phase!=-1) into a 2*Ncx × 2*Ncz fine mesh,
// then writes: count  num_cells_with_that_count
static void extractCellCountHistogram(const string &h5path, const string &tsvPath,
                                      int Ncx, int Ncz,
                                      double xmin, double xmax,
                                      double zmin, double zmax) {
  H5p::File file(h5path, "r");
  vector<float> px = file.read<vector<float>>("/Particles/x");
  vector<float> pz = file.read<vector<float>>("/Particles/z");
  vector<int8_t> pp = readInt8Dataset(h5path, "/Particles/phase");

  int nfx = 2 * Ncx, nfz = 2 * Ncz;
  double dx = (xmax - xmin) / nfx;
  double dz = (zmax - zmin) / nfz;
  vector<int> cellCount(nfx * nfz, 0);

  for (int i = 0; i < (int)px.size(); i++) {
    if (pp[i] == -1) continue;
    int ix = (int)((px[i] - xmin) / dx);
    int iz = (int)((pz[i] - zmin) / dz);
    if (ix < 0) ix = 0; if (ix >= nfx) ix = nfx - 1;
    if (iz < 0) iz = 0; if (iz >= nfz) iz = nfz - 1;
    cellCount[iz * nfx + ix]++;
  }

  // Build histogram: count -> how many cells have that count
  int maxCount = *max_element(cellCount.begin(), cellCount.end());
  vector<int> hist(maxCount + 1, 0);
  for (int c : cellCount) hist[c]++;

  ofstream out(tsvPath);
  out << "# particles_per_cell\tnum_cells\n";
  for (int c = 0; c <= maxCount; c++) {
    if (hist[c] > 0)
      out << c << "\t" << hist[c] << "\n";
  }
  out.close();
}

// Output 2D grid of per-fine-cell particle counts for heatmap plotting.
// Writes: x_centre  z_centre  count (one line per fine cell, blank line between rows)
static void extractCellCountGrid(const string &h5path, const string &tsvPath,
                                 int Ncx, int Ncz,
                                 double xmin, double xmax,
                                 double zmin, double zmax) {
  H5p::File file(h5path, "r");
  vector<float> px = file.read<vector<float>>("/Particles/x");
  vector<float> pz = file.read<vector<float>>("/Particles/z");
  vector<int8_t> pp = readInt8Dataset(h5path, "/Particles/phase");

  int nfx = 2 * Ncx, nfz = 2 * Ncz;
  double dx = (xmax - xmin) / nfx;
  double dz = (zmax - zmin) / nfz;
  vector<int> cellCount(nfx * nfz, 0);

  for (int i = 0; i < (int)px.size(); i++) {
    if (pp[i] == -1) continue;
    int ix = (int)((px[i] - xmin) / dx);
    int iz = (int)((pz[i] - zmin) / dz);
    if (ix < 0) ix = 0; if (ix >= nfx) ix = nfx - 1;
    if (iz < 0) iz = 0; if (iz >= nfz) iz = nfz - 1;
    cellCount[iz * nfx + ix]++;
  }

  ofstream out(tsvPath);
  for (int iz = 0; iz < nfz; iz++) {
    double zc = zmin + (iz + 0.5) * dz;
    for (int ix = 0; ix < nfx; ix++) {
      double xc = xmin + (ix + 0.5) * dx;
      out << xc << "\t" << zc << "\t" << cellCount[iz * nfx + ix] << "\n";
    }
    out << "\n";
  }
  out.close();
}


// Scan a directory for Particles*.h5, extract Nb_part from each
static void extractNbPartTimeseries(const string &dir, const string &tsvPath) {
  set<fs::path> files;
  for (auto &entry : fs::directory_iterator(dir)) {
    string fname = entry.path().filename().string();
    if (fname.find("Particles") == 0 && fname.find(".h5") != string::npos) {
      files.insert(entry.path());
    }
  }

  ofstream out(tsvPath);
  out << "# step\tNb_part\tNb_active\n";
  int step = 0;
  for (auto &p : files) {
    vector<int8_t> pp = readInt8Dataset(p.string(), "/Particles/phase");
    int nb = (int)pp.size();
    int nb_active = 0;
    for (auto ph : pp) if (ph != -1) nb_active++;

    // Extract step from filename: Particles00010.gzip.h5
    string fname = p.filename().string();
    size_t pos1 = fname.find("Particles") + 9;
    size_t pos2 = fname.find(".");
    if (pos1 != string::npos && pos2 != string::npos) {
      step = stoi(fname.substr(pos1, pos2 - pos1));
    }
    out << step << "\t" << nb << "\t" << nb_active << "\n";
  }
  out.close();
}

// ============================================================
//  Runner helpers
// ============================================================

static MdoodzSetup makeRotationSetup(bool multiphase = false) {
  MdoodzSetup setup = {
    .SetParticles = new SetParticles_ff{
      .SetHorizontalVelocity = RGSetHorizontalVelocity,
      .SetVerticalVelocity   = RGSetVerticalVelocity,
      .SetPhase              = multiphase ? RGSetPhaseMultiphase : RGSetPhase,
      .SetDensity            = RGSetDensity,
    },
    .SetBCs = new SetBCs_ff{
      .SetBCVx = SetPureShearBCVx,
      .SetBCVz = SetPureShearBCVz,
    },
  };
  return setup;
}

static void freeSetup(MdoodzSetup &s) {
  delete s.SetParticles;
  delete s.SetBCs;
}

static MdoodzSetup makeConvectionSetup() {
  MdoodzSetup setup = {
    .SetParticles = new SetParticles_ff{
      .SetHorizontalVelocity = CBSetHorizontalVelocity,
      .SetVerticalVelocity   = CBSetVerticalVelocity,
      .SetPhase              = RGSetPhase,
      .SetDensity            = RGSetDensity,
    },
    .SetBCs = new SetBCs_ff{
      .SetBCVx = SetPureShearBCVx,
      .SetBCVz = SetPureShearBCVz,
    },
  };
  return setup;
}

// ============================================================
//  Main entry points (called from main.cpp)
// ============================================================

void RunReseedingGuide() {
  // --- Mode comparison runs ---
  const char *modeConfigs[] = {
    "ReseedingGuide/CompareMode0.txt",
    "ReseedingGuide/CompareMode1.txt",
    "ReseedingGuide/CompareMode2.txt",
  };

  for (int m = 0; m < 3; m++) {
    MdoodzSetup setup = makeRotationSetup();
    RunMDOODZ((char*)modeConfigs[m], &setup);
    freeSetup(setup);
  }

  // --- Edge cases ---
  {
    MdoodzSetup setup = makeRotationSetup();
    RunMDOODZ((char*)"ReseedingGuide/EdgeLowPart.txt", &setup);
    freeSetup(setup);
  }
  {
    MdoodzSetup setup = makeRotationSetup();
    RunMDOODZ((char*)"ReseedingGuide/EdgeLowPart_Mode1.txt", &setup);
    freeSetup(setup);
  }
  {
    MdoodzSetup setup = makeRotationSetup();
    RunMDOODZ((char*)"ReseedingGuide/EdgeHighPart.txt", &setup);
    freeSetup(setup);
  }
  {
    MdoodzSetup setup = makeRotationSetup(true); // 3-phase
    RunMDOODZ((char*)"ReseedingGuide/EdgeMultiphase.txt", &setup);
    freeSetup(setup);
  }

  // --- High-density mode comparison (6×6 particles/cell) ---
  // Mode 0 omitted: CountPartCell_OLD has memory corruption with high particle density.
  const char *hdConfigs[] = {
    "ReseedingGuide/HighDensMode1.txt",
    "ReseedingGuide/HighDensMode2.txt",
  };
  for (int m = 0; m < 2; m++) {
    MdoodzSetup setup = makeRotationSetup();
    RunMDOODZ((char*)hdConfigs[m], &setup);
    freeSetup(setup);
  }

  // --- Closed-box convection test (no particles leave the domain) ---
  // Run each mode in a child process so exit(190) doesn't kill the runner
  const char *closedBoxConfigs[] = {
    "ReseedingGuide/ClosedBoxMode2.txt",
    "ReseedingGuide/ClosedBoxMode1.txt",
    "ReseedingGuide/ClosedBoxMode0.txt",
  };
  for (int m = 0; m < 3; m++) {
    printf("\n=== Running %s ===\n", closedBoxConfigs[m]);
    fflush(stdout);
    pid_t pid = fork();
    if (pid == 0) {
      MdoodzSetup setup = makeConvectionSetup();
      RunMDOODZ((char*)closedBoxConfigs[m], &setup);
      freeSetup(setup);
      _exit(0);
    } else if (pid > 0) {
      int status;
      waitpid(pid, &status, 0);
      if (WIFEXITED(status) && WEXITSTATUS(status) == 190)
        printf("=== %s CRASHED (exit 190) ===\n", closedBoxConfigs[m]);
    }
  }
}

void PlotReseedingGuide() {
  // --- Nb_part timeseries extraction ---
  const char *modes[] = {"Mode0", "Mode1", "Mode2"};
  for (int m = 0; m < 3; m++) {
    string dir = string("ReseedingGuide/") + modes[m];
    if (fs::exists(dir)) {
      extractNbPartTimeseries(dir, string("reseed_nbpart_") + modes[m] + ".dat");
    }
  }

  // --- Cell-level particle extraction for placement/deactivation plots ---
  // Target: a fine-mesh cell near the disk edge where reseeding happens
  // Fine-mesh dx = 0.5/20/2 = 0.0125; pick cell around (0.35, 0)
  double cellX = 0.35, cellHalf = 0.0125;
  double xlo = cellX - cellHalf, xhi = cellX + cellHalf;
  double zlo = -cellHalf, zhi = cellHalf;

  for (int m = 0; m < 3; m++) {
    string dir = string("ReseedingGuide/") + modes[m];
    // Initial state
    string h5_init = dir + "/Particles00000.gzip.h5";
    if (fs::exists(h5_init)) {
      extractParticlesForCell(h5_init, string("reseed_cell_init_") + modes[m] + ".dat",
                              xlo, xhi, zlo, zhi);
    }
    // After 25 steps
    string h5_mid = dir + "/Particles00025.gzip.h5";
    if (fs::exists(h5_mid)) {
      extractParticlesForCell(h5_mid, string("reseed_cell_mid_") + modes[m] + ".dat",
                              xlo, xhi, zlo, zhi);
      extractAllParticles(h5_mid, string("reseed_all_mid_") + modes[m] + ".dat");
    }
    // Final state
    string h5_final = dir + "/Particles00050.gzip.h5";
    if (fs::exists(h5_final)) {
      extractParticlesForCell(h5_final, string("reseed_cell_final_") + modes[m] + ".dat",
                              xlo, xhi, zlo, zhi);
      extractAllParticles(h5_final, string("reseed_all_final_") + modes[m] + ".dat");
    }
  }

  // --- Edge cases ---
  const char *edgeCases[] = {"EdgeLowPart", "EdgeLowPart_Mode1", "EdgeHighPart", "EdgeMultiphase"};

  // --- Extract all particle snapshots for Mode0/1/2 (for animated GIFs) ---
  const char *animModes[] = {"Mode0", "Mode1", "Mode2"};
  for (int am = 0; am < 3; am++) {
    string dir = string("ReseedingGuide/") + animModes[am];
    for (int step = 0; step <= 50; step += 1) {
      char h5name[128], datname[128], gridname[128];
      snprintf(h5name, sizeof(h5name), "%s/Particles%05d.gzip.h5", dir.c_str(), step);
      snprintf(datname, sizeof(datname), "reseed_anim_%s_%03d.dat", animModes[am], step);
      snprintf(gridname, sizeof(gridname), "reseed_grid_%s_%03d.dat", animModes[am], step);
      if (fs::exists(h5name)) {
        extractAllParticles(h5name, datname);
        extractCellCountGrid(h5name, gridname, 20, 20, -0.5, 0.5, -0.5, 0.5);
      }
    }
  }

  for (int e = 0; e < 4; e++) {
    string dir = string("ReseedingGuide/") + edgeCases[e];
    string h5_init = dir + "/Particles00000.gzip.h5";
    string h5_final = dir + "/Particles00050.gzip.h5";
    if (fs::exists(h5_init)) {
      extractAllParticles(h5_init, string("reseed_edge_init_") + edgeCases[e] + ".dat");
    }
    if (fs::exists(h5_final)) {
      extractAllParticles(h5_final, string("reseed_edge_final_") + edgeCases[e] + ".dat");
    }
    if (fs::exists(dir)) {
      extractNbPartTimeseries(dir, string("reseed_nbpart_") + edgeCases[e] + ".dat");
    }
  }

  // --- High-density mode comparison ---
  const char *hdModes[] = {"HighDensMode1", "HighDensMode2"};
  for (int m = 0; m < 2; m++) {
    string dir = string("ReseedingGuide/") + hdModes[m];
    if (fs::exists(dir)) {
      extractNbPartTimeseries(dir, string("reseed_nbpart_") + hdModes[m] + ".dat");
      string h5_final = dir + "/Particles00050.gzip.h5";
      if (fs::exists(h5_final)) {
        extractAllParticles(h5_final, string("reseed_hd_final_") + hdModes[m] + ".dat");
        extractCellCountHistogram(h5_final,
                                  string("reseed_cellhist_") + hdModes[m] + ".dat",
                                  20, 20, -0.5, 0.5, -0.5, 0.5);
      }
    }
  }
  // Also extract cell histograms for normal-density modes
  for (int m = 0; m < 3; m++) {
    string dir = string("ReseedingGuide/") + modes[m];
    string h5_final = dir + "/Particles00050.gzip.h5";
    if (fs::exists(h5_final)) {
      extractCellCountHistogram(h5_final,
                                string("reseed_cellhist_") + modes[m] + ".dat",
                                20, 20, -0.5, 0.5, -0.5, 0.5);
    }
  }

  // --- Closed-box convection: extract grids and Nb_part timeseries ---
  const char *cbModes[] = {"ClosedBoxMode0", "ClosedBoxMode1", "ClosedBoxMode2"};
  for (int m = 0; m < 3; m++) {
    string dir = string("ReseedingGuide/") + cbModes[m];
    if (fs::exists(dir)) {
      extractNbPartTimeseries(dir, string("reseed_nbpart_") + cbModes[m] + ".dat");
      for (int step = 0; step <= 3000; step += 30) {
        char h5name[128], gridname[128];
        snprintf(h5name, sizeof(h5name), "%s/Particles%05d.gzip.h5", dir.c_str(), step);
        snprintf(gridname, sizeof(gridname), "reseed_grid_%s_%03d.dat", cbModes[m], step);
        if (fs::exists(h5name)) {
          extractCellCountGrid(h5name, gridname, 10, 10, -0.5, 0.5, -0.5, 0.5);
        }
      }
    }
  }

  // --- Run gnuplot scripts ---
  const char *scripts[] = {
    "advection_primer.gnu",
    "reseed_placement.gnu",
    "deactivation_strategy.gnu",
    "nbpart_timeseries.gnu",
    "slot_usage.gnu",
    "edge_empty_cell.gnu",
    "edge_multiphase.gnu",
    "placement_zoom.gnu",
    "stress_comparison.gnu",
    "cell_uniformity.gnu",
    "advection_anim.gnu",
    "mode_comparison_anim.gnu",
    "closed_box_crash_anim.gnu",
  };
  for (auto *s : scripts) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "gnuplot -e \"outdir='../VISUAL_TESTS/img'\" %s", s);
    FILE *pipe = popen(cmd, "w");
    if (pipe) {
      fflush(pipe);
      pclose(pipe);
    }
  }
}
