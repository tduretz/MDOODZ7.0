/*
 * Extract temperature and velocity fields from MDOODZ HDF5 output
 * for gnuplot visualisation of the Blankenbach convection benchmark.
 *
 * Usage:
 *   ./extract_blankenbach BlankenBench/Output00500.gzip.h5
 *
 * Produces:
 *   BlankenBench/temperature.dat  — x(km) z(km) T/DeltaT  (pm3d format)
 *   BlankenBench/velocity.dat     — x(km) z(km) vx_norm vz_norm
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
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

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <Output.gzip.h5>\n", argv[0]);
    return 1;
  }

  const char *h5path = argv[1];

  // Derive output directory from input path
  std::string path(h5path);
  std::string outdir = ".";
  size_t sep = path.find_last_of("/\\");
  if (sep != std::string::npos) {
    outdir = path.substr(0, sep);
  }

  hid_t file = H5Fopen(h5path, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    fprintf(stderr, "ERROR: Cannot open %s\n", h5path);
    return 1;
  }

  // Read model parameters
  auto params = readDoubleDataset(file, "Model", "Params");
  int Nx  = (int)params[3];
  int Nz  = (int)params[4];
  int ncx = Nx - 1;
  int ncz = Nz - 1;

  // Read coordinates
  auto xc = readFloatDataset(file, "Model", "xc_coord");
  auto zc = readFloatDataset(file, "Model", "zc_coord");

  // Read fields
  auto T_flat  = readFloatDataset(file, "Centers", "T");
  auto Vx_flat = readFloatDataset(file, "VxNodes", "Vx");
  auto Vz_flat = readFloatDataset(file, "VzNodes", "Vz");

  H5Fclose(file);

  // Non-dimensionalise temperature: T_nd = (T - 273.15) / 1000
  const double T_top   = 273.15;
  const double DeltaT  = 1000.0;

  // Write temperature.dat (gnuplot pm3d: blank line between rows)
  std::string tpath = outdir + "/temperature.dat";
  FILE *fout = fopen(tpath.c_str(), "w");
  if (!fout) { fprintf(stderr, "Cannot write %s\n", tpath.c_str()); return 1; }
  fprintf(fout, "# x(km)  z(km)  T/DeltaT\n");
  for (int j = 0; j < ncz; j++) {
    for (int i = 0; i < ncx; i++) {
      double T_nd = (T_flat[j * ncx + i] - T_top) / DeltaT;
      fprintf(fout, "%.4f  %.4f  %.6f\n", xc[i] / 1e3f, zc[j] / 1e3f, T_nd);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
  printf("Wrote %s\n", tpath.c_str());

  // Interpolate Vx and Vz to cell centres
  // Vx: shape (Nz+1) x Nx  →  average columns to get ncx, skip first/last row for ncz
  // Vz: shape Nz x (Nx+1)  →  average rows to get ncz, skip first/last col for ncx
  std::vector<double> Vx_c(ncz * ncx), Vz_c(ncz * ncx);
  for (int j = 0; j < ncz; j++) {
    for (int i = 0; i < ncx; i++) {
      // Vx: row j+1 (skip ghost), cols i and i+1
      double vxL = Vx_flat[(j + 1) * Nx + i];
      double vxR = Vx_flat[(j + 1) * Nx + i + 1];
      Vx_c[j * ncx + i] = 0.5 * (vxL + vxR);
      // Vz: rows j and j+1, col i+1 (skip ghost)
      double vzB = Vz_flat[j * (Nx + 1) + i + 1];
      double vzT = Vz_flat[(j + 1) * (Nx + 1) + i + 1];
      Vz_c[j * ncx + i] = 0.5 * (vzB + vzT);
    }
  }

  // Find max speed for normalisation
  double maxSpeed = 0.0;
  for (int k = 0; k < ncz * ncx; k++) {
    double s = sqrt(Vx_c[k] * Vx_c[k] + Vz_c[k] * Vz_c[k]);
    if (s > maxSpeed) maxSpeed = s;
  }

  // Write velocity.dat
  std::string vpath = outdir + "/velocity.dat";
  fout = fopen(vpath.c_str(), "w");
  if (!fout) { fprintf(stderr, "Cannot write %s\n", vpath.c_str()); return 1; }
  fprintf(fout, "# x(km)  z(km)  vx_norm  vz_norm\n");
  for (int j = 0; j < ncz; j++) {
    for (int i = 0; i < ncx; i++) {
      int k = j * ncx + i;
      double vxn = (maxSpeed > 0) ? Vx_c[k] / maxSpeed : 0.0;
      double vzn = (maxSpeed > 0) ? Vz_c[k] / maxSpeed : 0.0;
      fprintf(fout, "%.4f  %.4f  %.6f  %.6f\n", xc[i] / 1e3f, zc[j] / 1e3f, vxn, vzn);
    }
  }
  fclose(fout);
  printf("Wrote %s\n", vpath.c_str());

  return 0;
}
