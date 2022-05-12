#include "math.h"
#include "mdoodz-private.h"
#include "stdio.h"


void PlotPhases(MdoodzInstance *instance, markers *particles) {
  const char *txtFileName = "phases.dat";

  const int Nx_part = particles->Nx_part;
  const int Nz_part = particles->Nz_part;
  const int Nb_part = particles->Nb_part;

  FILE *fp = fopen(txtFileName, "w");
  for (int np = 0; np < Nb_part; np += Nx_part * Nz_part) {
    const double x = particles->x[np];
    const double z = particles->z[np];
    const int phase = particles->phase[np];
    fprintf(fp, "%f\t%f\t%i\n", x, z, phase);
  }
  fclose(fp);

  FILE *GNUplotPipe = popen ("gnuplot -persistent", "w");
  fprintf(GNUplotPipe, "set view map'\n");
  fprintf(GNUplotPipe, "set dgrid3d\n");
  fprintf(GNUplotPipe, "set palette maxcolors %i\n", instance->model.Nb_phases);
  fprintf(GNUplotPipe, "set pm3d interpolate 4,4\n");
  fprintf(GNUplotPipe, "splot '%s' with pm3d\n", txtFileName);
  fflush(GNUplotPipe);
  exit(0);
}


void PlotResiduals(Nparams Nmodel, FILE *GNUplotPipe) {
  int   NumCommands       = 3;
  char *GNUplotCommands[] = {
          "set title \"Non-linear residuals\"",
          "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5",
          "set pointintervalbox 3"};
  for (int i = 0; i < NumCommands; i++) {
    fprintf(GNUplotPipe, "%s \n", GNUplotCommands[i]);
  }

  fprintf(GNUplotPipe, "plot '-' with linespoints linestyle 1\n");
  for (int i = 0; i < Nmodel.nit + 1; i++) {
    fprintf(GNUplotPipe, "%lf %lf \n", (double) i, log10(Nmodel.rx_rel[i]));//Write the data to a temporary file
  }
  fprintf(GNUplotPipe, "e\n");
  fflush(GNUplotPipe);
}