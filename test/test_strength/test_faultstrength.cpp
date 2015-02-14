/**
 * Test the strength accross the whole fault by outputting the
 * strength with depth and strike
 */
#include <iostream>
#include <algorithm>

#include "strength_profile.h"
#include "temperature_profile.h"
#include "creep_law.h"
#include "benzion_creep.h"
#include "side_creep.h"
#include "units.h"


/* Constants  */
const int nL(128), nD(32); // number of depth increments

const double arrhAmplitude( 6.31E-20 / ( pow(Units::Pa, 3) * Units::second) ), // Amplitude in the Arrhenius relation
  activationEnergy( 130*Units::kJ/Units::mole), // Activation energy at fault surface
  stressExponent( 3.0 ), // exponent to stress in arrhenius equation
  xBD( 7.5*Units::km ),  // creeping boundary width
  zBD( 10.5*Units::km ),
  Tsurface( ( 273.15 + 20 ) * Units::K ), // temperature of 283 K at surface
  dTdz( 20 * Units::K/Units::km ), // degrees K per km gradient
  cohesion( 60*Units::bar ) , // cohesion & max stress drop
  fs( 0.75 ) , // static coeff of friction
  dsigmadz( 180*(Units::bar/Units::km) ) , // rate of normal stress increase
  faultWidth( 10*Units::cm ), // width of the fault
  plateVelocity( 35*Units::mm/Units::year ), // velocity of fault loading
  cellLength( 70.0*Units::km/(double)nL ), // slip cell length
  cellHeight( 17.5*Units::km/(double)nD ); // slip cell height

//------------------------------------------------------------------------------
int main(int argc, char * argv[])
{

  std::cerr << "Program: " << argv[0] << " No. Args: " << argc << std::endl;

  /* Compute the loading strain rate on the fault */
  double loadingStrainRate = plateVelocity / faultWidth;

  /*Initialize the static strength profile */
  StaticStrengthProfile StatStrength( cohesion, fs, dsigmadz );

  /*Initialize the temperature profile */
  BackgroundTemperatureProfile T( Tsurface, dTdz );
  std::vector<double> faultBkgdT = T.faultgrid_co( nL, nD, cellHeight );

  /* Set up the creep parameters */
  std::vector<double> faultE(nL*nD, activationEnergy);
  std::vector<double> faultA(nL*nD, arrhAmplitude);
  std::vector<double> faultn(nL*nD, stressExponent);

  /* Apply a creep mask to the sides */
  sidecreepmask( faultA, faultE, faultn, faultBkgdT, nL, nD, cellLength, xBD );

  /* Calculate the creep strength */
  std::vector<double> creepStrength( nL*nD );
  for( unsigned int k=0; k<nL*nD; k++ ) {
    CreepLaw C0( faultA[k], faultE[k], faultn[k] );
    creepStrength[k] = C0.stress( loadingStrainRate, faultBkgdT[k] );
  }

  /* Calculate the creep mask used in BZ1996 */
  std::vector<double> cBZ96(nL*nD);
  faultcreep_bz96( cBZ96, nL, nD, cellLength, cellHeight, xBD, zBD,
                   loadingStrainRate, StatStrength( zBD ) );

  /* Calculate the BZ1996 creep strength */
  std::vector<double> creepStrengthBZ96( nL*nD );
  for( unsigned int k=0; k<nL*nD; k++ ) {
    CreepLaw C0( cBZ96[k], 0.0, 3.0 );
    creepStrengthBZ96[k] = C0.stress( loadingStrainRate, faultBkgdT[k] );
  }

  // /* Apply a creep mask to the sides */
  // sidecreepmask
  /** Change the output format */
  std::cout.setf(std::ios::scientific);
  std::cout.precision(15);

  /* Output the strength for the whole fault */
  for( int j=0; j<nD; j++ )
    for( int i=0; i<nL; i++ )
      {
        /* Calculate the positions */
	double xpos = ((double)i+0.5)*cellLength;
	double zpos = ((double)j+0.5)*cellHeight;

        /* Calculate the fault strength */
        double minStrength =  std::min( creepStrength[i*nD+j], StatStrength(zpos));
        double minStrengthBZ96 =  std::min( creepStrengthBZ96[i*nD+j], StatStrength(zpos));

        /* Calculate the strain rate at the fault strength */
        CreepLaw C0( faultA[i*nD+j], faultE[i*nD+j], faultn[i*nD+j] );
        double strainRate = C0.strainRate( minStrength, faultBkgdT[i*nD+j] );

        CreepLaw C1( cBZ96[i*nD+j], 0.0, 3.0 );
        double strainRateBZ96 = C1.strainRate( minStrengthBZ96, faultBkgdT[i*nD+j] );

	std::cout << xpos / Units::km // 1
		  << "\t " << zpos / Units::km // 2
		  << "\t " << faultA[i*nD+j] // 3
		  << "\t " << faultE[i*nD+j] // 4
                  << "\t " << faultBkgdT[i*nD+j]  // 5
		  << "\t " << StatStrength(zpos) / Units::MPa // 6
                  << "\t " << creepStrength[i*nD+j] / Units::MPa // 7
                  << "\t " << creepStrengthBZ96[i*nD+j] / Units::MPa // 8
                  << "\t " << minStrength / Units::MPa // 9
                  << "\t " << minStrengthBZ96 / Units::MPa // 10
                  << "\t " << strainRate / (1/Units::year) // 11
                  << "\t " << strainRateBZ96 / (1/Units::year) //12
		  << std::endl;
      }

}
/* End of test_faultstrength.cpp */
