#include "math.h"
#include "mdoodz-private.h"
#include "stdio.h"


void PlotPhases(MdoodzInstance *instance, markers *particles) {
  char *txtFileName = "phases.txt";
  const int cellParticles = particles->Nx_part * particles->Nz_part;
  const double resolution = 2;

  FILE *fp = fopen(txtFileName, "w");
  int particlesCount = 0;
  double xParticlesSum = 0.0;
  double zParticlesSum = 0.0;
  double phaseSum = 0.0;
  for (int np = 0; np < particles->Nb_part; np++) {
    particlesCount++;
    if (particlesCount < cellParticles * resolution) {
      const double x = particles->x[np];
      const double z = particles->z[np];
      const int phase = particles->phase[np];

      xParticlesSum += x;
      zParticlesSum += z;
      phaseSum += (double) phase;
    } else {
      const double x = xParticlesSum / (double) particlesCount;
      const double z = zParticlesSum / (double) particlesCount;
      const double phase = phaseSum / (double) particlesCount;
      fprintf(fp, "%f\t%f\t%f\n", x, z, phase);
      particlesCount = 0;
      xParticlesSum = 0.0;
      zParticlesSum = 0.0;
      phaseSum = 0.0;
      continue;
    }
  }
  fclose(fp);

  FILE *GNUplotPipe = popen ("gnuplot -persistent", "w");
  fprintf(GNUplotPipe, "load 'inferno.pal'\n");
  fprintf(GNUplotPipe, "set xrange[%f:%f]\n", instance->model.xmin, instance->model.xmax);
  fprintf(GNUplotPipe, "set yrange[%f:%f]\n", instance->model.zmin, instance->model.zmax);
  fprintf(GNUplotPipe, "plot '%s' with image\n", txtFileName);
  fflush(GNUplotPipe);
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