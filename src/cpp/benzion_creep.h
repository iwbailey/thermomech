#ifndef BENZION_CREEP_H_
#define BENZION_CREEP_H_
/**
 * Functions that define the distribution of activation energy over the fault
 *
 */

#include <vector>
#include <math.h>

/* Return the creep parameter c(x,z) from Ben-Zion (1996) */
double creep_bz96(
		  const double &x, // position from the fault edge
		  const double &vPl, // plate velocity or strain rate
		  const double &strengthz, // strength at the depth of the grid cell
		  const double &xBD ); // Width of the creep boundary

/* Return the creep coefficients c(x,z) for the entire fault, column ordered */
void faultcreep_bz96(
                     std::vector<double> &creepcoeffsBZ96, // OUT
                     const int &nL,
                     const int &nD,
                     const double &cellLength,
                     const double &cellHeight,
                     const double &xBD,
                     const double &zBD,
                     const double &vPl,
                     const double &strength_zBD );

#endif /*BENZION_CREEP_H_*/
