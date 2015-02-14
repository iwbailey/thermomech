#ifndef INPUT_PARAMETERS_H_
#define INPUT_PARAMETERS_H_

#include <vector>
#include <string.h>
#include "romberg_integration.h"
#include "fault_cooling.h"
#include "units.h"

/* Parameters for the integration */
int Romberg::JMAX = 20;
int Romberg::K = 5;
double Romberg::EPS = 1.0e-10;

/* Cooling model properties, must be defined upfront because static member

   - Density from PREM density of crustal rock kg/km^3
   - Width set arbitrarily, defines volume of fault heated by slip
   - Diffusivity can't remember origin
   - Specific heat capacity can't remember origin
*/
double FaultCoolingModel::density = 2.6*Units::g/(Units::cm*Units::cm*Units::cm);
double FaultCoolingModel::faultWidth = 10.0*Units::cm;
double FaultCoolingModel::diffusivity = 1.E-2*Units::cm*Units::cm/Units::second;
//double FaultCoolingModel::diffusivity = 1.E-3*Units::cm*Units::cm/Units::second;
double FaultCoolingModel::specHeat = 790*Units::J/(Units::kg*Units::K);

//----------------------------------------
/**  Set parameters used in main function
*/
namespace In
{
  /* Fault properties */
  const double
    cellLength( 70.0*Units::km/128.0 ), // slip cell length
    cellHeight( 17.5*Units::km/32.0 ), // slip cell height
    cellDepth( 10.0*Units::km ), // Depth of the center of the slip cell
    faultWidth( 10*Units::cm ), // width of the fault for strain calc
    rigidity( 30*Units::GPa ); // rigidity of the halfspace

  /* Loading and creep parameters */
  const double plateVelocity( 35*Units::mm/Units::year ), // plate velocity
    loadingStrainRate( plateVelocity/faultWidth ) , // strain rate of fault loading
    arrhAmplitude( 6.31E-20 / ( pow(Units::Pa, 3) * Units::second) ), // Amplitude in the Arrhenius relation
    activationEnergy( 130*Units::kJ/Units::mole), // Activation energy at fault surface
    stressExponent( 3.0 ); // exponent to stress in arrhenius equation

  /* Background temperature profile */
  const double Tsurface( 273.15*Units::K ),  // Temperature at surface
    Tgradient( 20.0 * Units::K/Units::km ); // Change in temperature with depth

  /* Strength and initial stress parameters */
  const double
    tau0( 6.0 * Units::MPa ), // cohesion
    fs(0.75), // coefficient of friction
    dsigmadz( 18.0 * Units::MPa/Units::km ), // effective normal stress gradient
    dynStrengthDrop( 4.0*Units::MPa ), // tau_static - tau_dynamic
    dosCoef(1.25);

  /* Time taken for 1 m of seismic slip */
  const double slipVelocity(6.0 * Units::km/Units::second);

  /* Max amount of creep allowed in one time step */
  const double maxCreepSlip( 10.0*Units::cm );

  /* Algorithm time constraints */
  const double maxTimeStep = 1.0*Units::day;
  const double minTimeStep = 1.0*Units::minute;
  const int nTimeMax = 1e6; // maximum number of time steps
  const double maxTime = 10*Units::year;

}

#endif /* PARAMETERS_H_ */
