#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <vector>
#include <string.h>
#include <math.h>
#include "units.h"


//----------------------------------------
/**  Set parameters used in main function
*/
namespace In
{
  /* Fault properties */
  const int nL(128), nD(32); // number of cells along strike and down-dip

  const double cellLength( 70.0*Units::km/(double)nL ), // slip cell length
    cellHeight( 17.5*Units::km/(double)nD ), // slip cell height
    faultWidth( 10*Units::cm ), // width of the fault for strain calc
    rigidity( 30*Units::GPa ); // rigidity of the halfspace

  /* Loading and creep parameters */
  const double plateVelocity( 35*Units::mm/Units::year ), // plate velocity
    loadingStrainRate( plateVelocity/faultWidth ) , // strain rate of fault loading
    arrhAmplitude( 6.31E-20 / ( pow(Units::Pa, 3) * Units::second) ), // Amplitude in the Arrhenius relation
    activationEnergy( 130*Units::kJ/Units::mole), // Activation energy at fault surface
    stressExponent( 3.0 ), // exponent to stress in arrhenius equation
    xBD( 7.5*Units::km ); // creeping boundary width

  /* Background temperature profile */
  const double Tsurface( 273.15*Units::K ),  // Temperature at surface
    Tgradient( 20.0 * Units::K/Units::km ); // Change in temperature with depth

  /* Strength and initial stress parameters */
  const double
    tau0( 6.0 * Units::MPa ), // cohesion
    fs(0.75), // coefficient of friction
    dsigmadz( 18.0 * Units::MPa/Units::km ), // effective normal stress gradient
    dynStrengthDrop( 1.0*Units::MPa ), // tau_static - tau_dynamic
    dosCoef(1.25);

  /* Distribution of tau_s - tau_d, should be in Pa */
  //const std::string ifile_strengthdrops("../../inputs/stressdrops_rwalk.txt");
  const std::string ifile_strengthdrops("../../inputs/stressdrops_unif.txt");
  // const std::string ifile_strengthdrops("../../inputs/stressdrops_frac.txt");

  /* Time taken for 1 m of slip */
  const double slipVelocity(6.0 * Units::km/Units::second);

  /* Algorithm time constraints */
  const double maxTimeStep = 1.0*Units::day;
  const double minTimeStep = 1.0*Units::minute;
  const int nTimeMax = 1e6; // maximum number of time steps
  const double maxTime = 100.0*Units::year;

  /* Output file names */
  const std::string ofilesuffix_initStress = "_initstress.txt";
  const std::string ofilesuffix_finalStress = "_finalstress.txt";
  const std::string ofilesuffix_creepSlip = "_creepslip.txt";
  const std::string ofilesuffix_eqSlip = "_eqkslip.txt";
  const std::string ofilesuffix_slipdef = "_finalslipdef.txt";

}

#endif /* PARAMETERS_H_ */
