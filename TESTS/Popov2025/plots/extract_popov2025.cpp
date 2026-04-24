// Extract per-step centre-grid fields from a Popov2025 test HDF5 output into
// plain .dat files suitable for gnuplot. Follows the same pattern as
// TESTS/BlankenBench/extract_fields.cpp but iterates over multiple step files
// and emits one .dat per (field, step) plus a consolidated per-field .dat for
// centre-cell trajectories.
//
// Usage:
//   ./extract_popov2025 <test_output_dir> [--trajectory]
//
// <test_output_dir> is the writer_subfolder of the test (e.g.
// "Popov0D_VolumetricExtension"). The tool scans for Output*.gzip.h5 files in
// that directory and writes:
//
//   <test_output_dir>/trajectory.dat   (columns: step t sII P eII_pl divu_pl)
//   <test_output_dir>/sII_stepXXXXX.dat
//   <test_output_dir>/P_stepXXXXX.dat
//   ...
//
// For the 2D localisation tests, the per-step .dat files are the pm3d-format
// field maps consumed by TESTS/Popov2025/plots/fig6_regularization.gp.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <dirent.h>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <algorithm>
#include <hdf5.h>

static std::vector<float> readFloatDataset(hid_t file, const char *group, const char *name) {
  hid_t grp   = H5Gopen(file, group, H5P_DEFAULT);
  hid_t dset  = H5Dopen(grp, name, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  std::vector<float> buf(dims[0]);
  H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(grp);
  return buf;
}

static std::vector<double> readDoubleDataset(hid_t file, const char *group, const char *name) {
  hid_t grp   = H5Gopen(file, group, H5P_DEFAULT);
  hid_t dset  = H5Dopen(grp, name, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  std::vector<double> buf(dims[0]);
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(grp);
  return buf;
}

static std::vector<std::string> listStepFiles(const std::string &dir) {
  std::vector<std::string> out;
  DIR *d = opendir(dir.c_str());
  if (!d) return out;
  struct dirent *e;
  while ((e = readdir(d)) != NULL) {
    const std::string name = e->d_name;
    if (name.rfind("Output", 0) == 0 && name.find(".gzip.h5") != std::string::npos) {
      out.push_back(dir + "/" + name);
    }
  }
  closedir(d);
  std::sort(out.begin(), out.end());
  return out;
}

static void writePm3d(const std::string &path,
                      const std::vector<float> &xc, const std::vector<float> &zc,
                      const std::vector<float> &field, int ncx, int ncz,
                      const char *header) {
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp) { fprintf(stderr, "Cannot write %s\n", path.c_str()); return; }
  fprintf(fp, "# %s\n", header);
  for (int j = 0; j < ncz; ++j) {
    for (int i = 0; i < ncx; ++i) {
      fprintf(fp, "%.6g  %.6g  %.6g\n", (double)xc[i], (double)zc[j],
              (double)field[j * ncx + i]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

static double centreValue(const std::vector<float> &field, int ncx, int ncz) {
  const int cx = (ncx - 1) / 2;
  const int cz = (ncz - 1) / 2;
  return field[cz * ncx + cx];
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <test_output_dir>\n", argv[0]);
    return 1;
  }
  const std::string dir = argv[1];

  const auto files = listStepFiles(dir);
  if (files.empty()) {
    fprintf(stderr, "No Output*.gzip.h5 files in %s\n", dir.c_str());
    return 1;
  }

  const std::string trajPath = dir + "/trajectory.dat";
  FILE *traj = fopen(trajPath.c_str(), "w");
  if (!traj) { fprintf(stderr, "Cannot write %s\n", trajPath.c_str()); return 1; }
  fprintf(traj, "# step  time[s]  sII[Pa]  P[Pa]  eII_pl[1/s]  divu_pl[1/s]\n");

  for (size_t fi = 0; fi < files.size(); ++fi) {
    hid_t fd = H5Fopen(files[fi].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fd < 0) { fprintf(stderr, "Cannot open %s\n", files[fi].c_str()); continue; }

    auto params = readDoubleDataset(fd, "Model", "Params");
    const int Nx  = (int)params[3];
    const int Nz  = (int)params[4];
    const int ncx = Nx - 1;
    const int ncz = Nz - 1;
    const double t = params[0];
    // Parse step from filename "Output<NNNNN>.gzip.h5" (params[1] is Lx,
    // not the step counter as originally assumed).
    int step = 0;
    {
      const size_t p1 = files[fi].find("Output");
      if (p1 != std::string::npos) {
        step = std::atoi(files[fi].c_str() + p1 + 6);
      }
    }

    auto xc = readFloatDataset(fd, "Model", "xc_coord");
    auto zc = readFloatDataset(fd, "Model", "zc_coord");

    // MDOODZ doesn't write sII directly — reconstruct from Centers/sxxd,
    // Centers/szzd, and Vertices/sxz (averaged from four surrounding
    // vertices).
    auto sxxd = readFloatDataset(fd, "Centers", "sxxd");
    auto szzd = readFloatDataset(fd, "Centers", "szzd");
    auto sxz  = readFloatDataset(fd, "Vertices", "sxz");
    std::vector<float> sII(ncx * ncz, 0.0f);
    for (int j = 0; j < ncz; ++j) {
      for (int i = 0; i < ncx; ++i) {
        const int idxC = j * ncx + i;
        const double sxz_c = 0.25 * ((double)sxz[j * Nx + i]
                                   + (double)sxz[j * Nx + (i + 1)]
                                   + (double)sxz[(j + 1) * Nx + i]
                                   + (double)sxz[(j + 1) * Nx + (i + 1)]);
        const double inv2 = 0.5 * ((double)sxxd[idxC] * (double)sxxd[idxC]
                                 + (double)szzd[idxC] * (double)szzd[idxC])
                          + sxz_c * sxz_c;
        sII[idxC] = (float)std::sqrt(inv2 > 0.0 ? inv2 : 0.0);
      }
    }
    auto P      = readFloatDataset(fd, "Centers", "P");
    auto eIIpl  = readFloatDataset(fd, "Centers", "eII_pl");
    auto divupl = readFloatDataset(fd, "Centers", "divu_pl");

    fprintf(traj, "%d  %.6e  %.6e  %.6e  %.6e  %.6e\n",
            step, t,
            (double)centreValue(sII, ncx, ncz),
            (double)centreValue(P,   ncx, ncz),
            (double)centreValue(eIIpl,  ncx, ncz),
            (double)centreValue(divupl, ncx, ncz));

    // Per-step pm3d maps (only useful for the 2D tests).
    char stepTag[32];
    std::snprintf(stepTag, sizeof(stepTag), "_step%05d", step);

    writePm3d(dir + "/sII"    + stepTag + ".dat", xc, zc, sII,    ncx, ncz, "x z sII[Pa]");
    writePm3d(dir + "/P"      + stepTag + ".dat", xc, zc, P,      ncx, ncz, "x z P[Pa]");
    writePm3d(dir + "/eII_pl" + stepTag + ".dat", xc, zc, eIIpl,  ncx, ncz, "x z eII_pl[1/s]");
    writePm3d(dir + "/divu_pl"+ stepTag + ".dat", xc, zc, divupl, ncx, ncz, "x z divu_pl[1/s]");

    H5Fclose(fd);
  }

  fclose(traj);
  printf("Wrote %s (and per-step pm3d maps) for %zu output files.\n",
         trajPath.c_str(), files.size());
  return 0;
}
