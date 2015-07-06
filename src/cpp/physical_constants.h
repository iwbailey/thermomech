#ifndef PHYSICAL_CONSTANTS_H_
#define PHYSICAL_CONSTANTS_H_
/**
   Namespace defining physical constants necessary for the equations used
 */
#include "units.h"

namespace PhysicalConstants
{
  /* universal gas const. J/(mol*K) */
  static const double R_g = 8.3144 * Units::J / (Units::mole * Units::K );
}

#endif /*PHYSICAL_CONSTANTS_H_*/
