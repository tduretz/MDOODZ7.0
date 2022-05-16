extern "C" {
#include "visual-tests.h"
}

int main() {
  MdoodzInstance shearTemplateStyle1 = CreateShearTemplateInstance(1);
  MdoodzInstance shearTemplate       = CreateShearTemplateInstance(0);
  MdoodzInstance riftingPauline      = CreateRiftingPaulineInstance();
  MdoodzInstance topoBenchCase       = CreateTopoBenchCase1Instance();

  RunMDOODZ(&shearTemplateStyle1);
  RunMDOODZ(&shearTemplate);
  RunMDOODZ(&riftingPauline);
  RunMDOODZ(&topoBenchCase);
}