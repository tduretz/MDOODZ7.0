/*
 * Standalone runner for the closed-box convection reseeding test.
 * Runs Mode 2, Mode 1, Mode 0 in separate child processes (via fork)
 * so that exit(190) from Nb_part_max overflow doesn't kill the runner.
 * Then extracts data and generates the GIF.
 *
 * Build: linked against the same libraries as visualtests
 * Usage: ./closed_box_runner
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

extern "C" {
#include "mdoodz.h"
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace fs = std::filesystem;
using namespace std;

// ============================================================
//  Velocity callbacks — closed-box convection roll
// ============================================================

static const double CONV_AMP = 200.0;

static int CBSetPhase(MdoodzInput *input, Coordinates coords) {
  (void)input; (void)coords;
  return 0;
}

static double CBSetDensity(MdoodzInput *input, Coordinates coords, int phase) {
  return input->materials.rho[phase];
}

static double CBSetHVel(MdoodzInput *input, Coordinates coords) {
  double x = coords.x * input->scaling.L;
  double z = coords.z * input->scaling.L;
  return CONV_AMP * sin(M_PI * (x + 0.5)) * cos(M_PI * (z + 0.5)) / input->scaling.V;
}

static double CBSetVVel(MdoodzInput *input, Coordinates coords) {
  double x = coords.x * input->scaling.L;
  double z = coords.z * input->scaling.L;
  return -CONV_AMP * cos(M_PI * (x + 0.5)) * sin(M_PI * (z + 0.5)) / input->scaling.V;
}

static MdoodzSetup makeSetup() {
  MdoodzSetup setup = {
    .SetParticles = new SetParticles_ff{
      .SetHorizontalVelocity = CBSetHVel,
      .SetVerticalVelocity   = CBSetVVel,
      .SetPhase              = CBSetPhase,
      .SetDensity            = CBSetDensity,
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

// ============================================================
//  Helper: read int8 datasets
// ============================================================

static vector<int8_t> readInt8(const string &h5, const string &ds) {
  H5::H5File f(h5, H5F_ACC_RDONLY);
  H5::DataSet d = f.openDataSet(ds);
  size_t n = d.getSpace().getSelectNpoints();
  vector<int8_t> v(n);
  d.read(v.data(), H5::PredType::NATIVE_INT8);
  return v;
}

// ============================================================
//  Extractors
// ============================================================

static void extractNbPart(const string &dir, const string &out) {
  set<fs::path> files;
  for (auto &e : fs::directory_iterator(dir)) {
    string f = e.path().filename().string();
    if (f.find("Particles") == 0 && f.find(".h5") != string::npos)
      files.insert(e.path());
  }
  ofstream o(out);
  o << "# step\tNb_part\tNb_active\n";
  for (auto &p : files) {
    auto pp = readInt8(p.string(), "/Particles/phase");
    int nb = (int)pp.size(), act = 0;
    for (auto ph : pp) if (ph != -1) act++;
    string fn = p.filename().string();
    int step = stoi(fn.substr(9, fn.find(".") - 9));
    o << step << "\t" << nb << "\t" << act << "\n";
  }
}

static void extractGrid(const string &h5, const string &out,
                         int Ncx, int Ncz,
                         double xmin, double xmax, double zmin, double zmax) {
  H5p::File file(h5, "r");
  auto px = file.read<vector<float>>("/Particles/x");
  auto pz = file.read<vector<float>>("/Particles/z");
  auto pp = readInt8(h5, "/Particles/phase");

  int nfx = 2*Ncx, nfz = 2*Ncz;
  double dx = (xmax-xmin)/nfx, dz = (zmax-zmin)/nfz;
  vector<int> cc(nfx*nfz, 0);
  for (int i = 0; i < (int)px.size(); i++) {
    if (pp[i] == -1) continue;
    int ix = (int)((px[i]-xmin)/dx), iz = (int)((pz[i]-zmin)/dz);
    if (ix<0) ix=0; if (ix>=nfx) ix=nfx-1;
    if (iz<0) iz=0; if (iz>=nfz) iz=nfz-1;
    cc[iz*nfx+ix]++;
  }
  ofstream o(out);
  for (int iz=0; iz<nfz; iz++) {
    double zc = zmin + (iz+0.5)*dz;
    for (int ix=0; ix<nfx; ix++) {
      o << xmin+(ix+0.5)*dx << "\t" << zc << "\t" << cc[iz*nfx+ix] << "\n";
    }
    o << "\n";
  }
}

// ============================================================
//  Main
// ============================================================

int main() {
  // Run simulations: Mode 2 first (safe), then 1, then 0
  // Each runs in a child process so that exit(190) doesn't kill us
  const char *configs[] = {
    "ReseedingGuide/ClosedBoxMode2.txt",
    "ReseedingGuide/ClosedBoxMode1.txt",
    "ReseedingGuide/ClosedBoxMode0.txt",
  };
  for (int i = 0; i < 3; i++) {
    printf("\n=== Running %s ===\n", configs[i]);
    fflush(stdout);

    pid_t pid = fork();
    if (pid == 0) {
      // Child: run the simulation
      MdoodzSetup s = makeSetup();
      RunMDOODZ((char*)configs[i], &s);
      freeSetup(s);
      _exit(0);
    } else if (pid > 0) {
      // Parent: wait for child
      int status;
      waitpid(pid, &status, 0);
      if (WIFEXITED(status)) {
        int code = WEXITSTATUS(status);
        if (code == 190 || code == 1) {
          printf("=== %s CRASHED (exit %d — Nb_part_max exceeded) ===\n", configs[i], code);
        } else if (code != 0) {
          printf("=== %s exited with code %d ===\n", configs[i], code);
        }
      } else if (WIFSIGNALED(status)) {
        printf("=== %s killed by signal %d ===\n", configs[i], WTERMSIG(status));
      }
    } else {
      perror("fork failed");
      return 1;
    }
  }

  // Extract data — grid is 10x10 (Nx=Nz=11)
  const int NCX = 10, NCZ = 10;
  const int NT = 3000, WSTEP = 30;
  const char *cbm[] = {"ClosedBoxMode0", "ClosedBoxMode1", "ClosedBoxMode2"};
  for (int m = 0; m < 3; m++) {
    string dir = string("ReseedingGuide/") + cbm[m];
    if (!fs::exists(dir)) continue;
    extractNbPart(dir, string("reseed_nbpart_") + cbm[m] + ".dat");
    for (int step = 0; step <= NT; step += WSTEP) {
      char h5[128], grd[128];
      snprintf(h5, sizeof(h5), "%s/Particles%05d.gzip.h5", dir.c_str(), step);
      snprintf(grd, sizeof(grd), "reseed_grid_%s_%03d.dat", cbm[m], step);
      if (fs::exists(h5)) extractGrid(h5, grd, NCX, NCZ, -0.5, 0.5, -0.5, 0.5);
    }
  }

  // Generate GIF
  system("gnuplot -e \"outdir='../VISUAL_TESTS/img'\" closed_box_crash_anim.gnu");

  printf("\nDone! GIF at VISUAL_TESTS/img/closed_box_crash_anim.gif\n");
  return 0;
}
