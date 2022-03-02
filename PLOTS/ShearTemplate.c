#include "header_MDOODZ.h"
#include "stdio.h"
#define FILENAME "Output00010.gzip.h5"

int main() {
  char setupFileName[] = "ShearTemplatePlot.txt";
  const char *args[] = {
      "Some string",
      setupFileName,
  };
  RunMDOODZ(2, (char **)args);

  rename(FILENAME, "ShearTemplateResult.gzip.h5");
}