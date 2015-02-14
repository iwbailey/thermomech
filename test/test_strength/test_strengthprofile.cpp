/**
 * Test the strength profile and creep law classes by outputting the
 * strength with depth
 */
#include <iostream>
#include <algorithm>

#include "strength_profile.h"
#include "temperature_profile.h"
#include "creep_law.h"
#include "benzion_creep.h"
#include "units.h"

/* Constants  */
const double
arrhAmplitude( 6.31E-20 / ( pow(Units::Pa, 3) * Units::second) ), // Amplitude in the Arrhenius relation
  activationEnergy( 130*Units::kJ/Units::mole), // Activation energy at fault surface
  Tsurface( ( 273.15 + 20 ) * Units::K ), // temperature of 283 K at surface
  dTdz( 20 * Units::K/Units::km ), // degrees K per km gradient
  cohesion( 60*Units::bar ) , // cohesion & max stress drop
  fs( 0.75 ) , // static coeff of friction
  sigmaNgrad( 280*(Units::bar/Units::km) ) , // rate of normal stress increase
  porePgrad( 100*(Units::bar/Units::km) ), // rate of pore pressure increase
  faultWidth( 10*Units::cm ), // width of the fault
  plateVelocity( 35*Units::mm/Units::year ), // velocity of fault loading
  z0(0.0), // minimum depth for depth profile
  dz((17.5*Units::km)/32), // increments for depth
  zBD_bz96(10.0*Units::km); // Brittle-Ductile transition depth for BZ 1996

//const int nz(70); // number of depth increments
const int nz(32); // number of depth increments


//------------------------------------------------------------------------------
int main(int argc, char * argv[]) {

  std::cerr << "Program: " << argv[0] << " No. Args: " << argc << std::endl;

  /*Initialize the strength profile */
  StaticStrengthProfile Strength( cohesion, fs, sigmaNgrad, porePgrad );

  /*Initialize the creep law */
  CreepLaw C( arrhAmplitude, activationEnergy );

  /*Initialize the temperature profile */
  BackgroundTemperatureProfile T( Tsurface, dTdz );

  /*Compute the strain rate */
  double strainrate = plateVelocity/faultWidth;

  /*Compute fault depth*/
  double faultDepth = (double)(nz+1)*dz;

  /** Change the output format */
  std::cout.setf(std::ios::scientific);
  std::cout.precision(15);

  /*Loop through all depths */
  for( int i=0; i<nz; i++ ){

    double depth = (double)i*dz;

    /* Get temperature*/
    double temperature = T(depth);

    /* Get strength */
    double strengthStatic = Strength(depth);

    /* Get strength in BZ 1996 creep model */
    double cBZ1996 = creep_bz96( faultDepth-depth, strainrate,
                                 Strength(zBD_bz96), faultDepth-zBD_bz96 );
    CreepLaw Cbz96( cBZ1996, 0.0, 3);

    /* Get stress required for creeep at plate velocity */
    double strengthCreep = C.stress( strainrate, temperature );
    std::cout << depth / Units::km
	      << "\t " << strengthStatic / Units::MPa
	      << "\t " << strengthCreep /Units::MPa
	      << "\t " << std::min( strengthStatic, strengthCreep ) /Units::MPa
              << "\t " << Cbz96.stress( strainrate, temperature ) / Units::MPa
	      << std::endl;
  }
}
//======================================================================
//######################################################################
//######################################################################
