#ifndef FAULT_COOLING_H_
#define FAULT_COOLING_H_

#include <math.h>
#include "romberg_integration.h"

//------------------------------------------------------------------------------
/**
    A halfspace with constant heat decay properties surrounds a fault
    with constant finite width.

    Cardwell et al (1978)
 */
class FaultCoolingModel{

 public:
  /* All members are static and must be defined at in the initial part
     of the main source file.  Use static because should have only one
     decay function, and makes passing a function for integration
     easier. */
  static double faultWidth; // default 10 cm
  static double density; // default is 2.6 g/cm^3
  static double diffusivity; // default is 0.01 cm^2/s
  static double specHeat; // default is 790 J/(kg K)

  /* Temperature generated immediately from heat rate spike*/
  double temperature(
		     const double &heatRate, // heat rate in J/s
		     const double &area ); // area of heat generation in m^2

  /* Temperature change at time t from constant heat rate between t1 and t2*/
  double temperature(
		     const double &heatRate, // heat rate in J/s
		     const double &area, // area of heat generation
		     const double &t1, // start time for heat rate
		     const double &t2, // end time for heat rate
		     const double &tNow ); // time to calculate temperature at

  /* Temperature decay function for generated heat  */
  double integTempDecay(
			const double &tSlipStart, // start time of heat generation
			const double &tSlipEnd, // end time of heat generation
			const double &tNow ); // measurement time

private:
  /* Romberg object for performing the integration */
  Romberg R;

  /* Temperature decay function for a unit pulse of heat, at
     t=timeSinceSlip ago */
  static double tempDecay( double tSinceSlip );

};
//------------------------------------------------------------------------------
#endif /*FAULT_COOLING_H_*/
