/**
 * test_loading_noheat.cpp ---
 *
 *  Test the stress build up on a fault that has creep but no earthquakes and no
 *  heat generation
 *
 *
 * Filename: test_loading_noheat.cpp
 * Description:
 * Author: Iain W. Bailey
 * Created: Wed Sep 25 11:31:04 2013 (-0400)
 * Version: 1
 * Last-Updated: Sun May  4 20:51:14 2014 (-0400)
 *           By: Iain W. Bailey
 *     Update #: 63
 *

 * Change Log:
 *
 *
 *
 */
#include <vector>
#include <iostream>

#include "units.h"
#include "temperature_profile.h"
#include "strength_profile.h"
#include "side_creep.h"
#include "creep_law.h"
#include "stiffness_matrix.h"
#include "fault.h"

//------------------------------------------------------------------------------
int main(int argc, char * argv[])
{

  std::cerr << "Program: " << argv[0] << " No. Args: " << argc << std::endl;

  //----------------------------------------
  /**  Set parameters */
  /* Fault properties */
  const int
    nL(128), // number of cells along-strike
    nD(32); // number of cells down-dip

  const double
    cellLength( 70.0*Units::km/(double)nL ), // slip cell length
    cellHeight( 17.5*Units::km/(double)nD ), // slip cell height
    faultWidth( 10*Units::cm ), // width of the fault for strain calc
    rigidity( 30*Units::GPa ); // rigidity of the halfspace

  /* Loading and creep parameters */
  const double
    plateVelocity( 35*Units::mm/Units::year ), // plate velocity
    loadingStrainRate( plateVelocity/faultWidth ) , // strain rate of fault loading
    arrhAmplitude( 6.31E-20 / ( pow(Units::Pa, 3) * Units::second) ), // Amplitude in the Arrhenius relation
    activationEnergy( 130*Units::kJ/Units::mole), // Activation energy at fault surface
    stressExponent( 3.0 ), // exponent to stress in arrhenius equation
    xDB( 7.5*Units::km ); // creeping boundary width

  /* Strength and initial stress parameters */
  const double
    tau0( 6.0 * Units::MPa ), // cohesion
    fs(0.75), // coefficient of friction
    dsigmadz( 18.0 * Units::MPa/Units::km ), // effective normal stress gradient
    dynStrengthDrop( 1.0*Units::MPa ), // tau_static - tau_dynamic
    dosCoef(1.25);

  const double timeStep = 1.0*Units::year;
  const int nTime = 100; // number of time steps to load

  /* Background temperature */
  BackgroundTemperatureProfile T( 273.15*Units::K, 20.0 * Units::K/Units::km );

  /* Static strength profile */
  StaticStrengthProfile StatStrength( tau0, fs, dsigmadz );

  //----------------------------------------

  /* Compute the stiffness matrix */
  StiffnessMatrix K( cellLength, nL, cellHeight, nD, rigidity );

  // /* Compute the base strength profile */
  // CreepLaw C0( arrhAmplitude, activationEnergy, stressExponent );
  // std::vector<double> tempProfile(nD);
  // double peakStrength = 0.0;
  // for( int j=0; j<nD; j++ )
  //   {
  //     double depth = ((double)j+0.5)*cellHeight;
  //     double strengthDyn = C0.stress( loadingStrainRate, T(depth) );
  //     double minStrength = std::min( StatStrength(depth), strengthDyn );
  //     peakStrength = std::max( peakStrength, minStrength);
  //     tempProfile[j] = T(depth);
  //   }

  /* Get the  background temperature */
  std::vector<double> faultBkgdT = T.faultgrid_co( nL, nD, cellHeight );

  /* Set up the creep parameters over entire fault */
  std::vector<double> faultE(nL*nD, activationEnergy);
  std::vector<double> faultA(nL*nD, arrhAmplitude);
  std::vector<double> faultn(nL*nD, stressExponent );

  /* Set creep mask */
  sidecreepmask( faultA, faultE, faultn, faultBkgdT, nL, nD, cellLength, xDB );

  /*Initialize the background temperature and strength over entire fault */
  std::vector<double> initStress(nL*nD), taus(nL*nD), taud(nL*nD);
  for( int j=0; j<nD; j++ )
    {
      double depth = ((double)j + 0.5)*cellHeight;
      for( int i=0; i<nL; i++ )
	{
	  CreepLaw C( faultA[i*nD+j], faultE[i*nD+j], faultn[i*nD+j] );
	  initStress[i*nD+j] = std::min( StatStrength(depth),
					 C.stress( loadingStrainRate,
                                                   faultBkgdT[i*nD+j] ));
	  taus[i*nD+j] = StatStrength(depth);
	  taud[i*nD+j] = taus[i*nD+j] - dynStrengthDrop;
	}
    }


  /** Initialize the fault */
  Fault F( nL, nD, cellLength, cellHeight, faultWidth, taus, taud, dosCoef,
	   faultA, faultn, faultE, K, initStress );

  /** Load the fault */
  std::vector<double> stress(nL*nD), creepVel(nL*nD);
  for( int i=0; i<nTime; i++ )
    {
      std::cerr << "time step " << i << std::endl;

      /* Get current state of stress and temperature */
      F.getStress( stress );

      /* Get the creep velocity on the fault */
      F.getCreepVelocity( creepVel, stress, faultBkgdT);

      /* Apply the loading */
      F.loadFault( plateVelocity, creepVel, timeStep );
    }

  /** Calculate stress */
  F.getStress( stress );
  F.getCreepVelocity( creepVel, stress, faultBkgdT);

  /* Output stress and slip deficit to terminal */
  std::cout.setf(std::ios::scientific);
  std::cout.precision(7);

  for( int j=0; j<nD; j++ )
    for( int i=0; i<nL; i++ )
      std::cout << ((double)i + 0.5)*cellLength / Units::km
		<< "\t " << ((double)j + 0.5)*cellHeight / Units::km
		<< "\t " << F.getCellStress( i, j) / Units::MPa
		<< "\t " << F.getCellSlipDeficit( i, j) / Units::m
                << "\t " << faultBkgdT[i*nD+j]
                << "\t " << creepVel[i*nD+j] / (Units::mm / Units::year )
                << "\t " << faultA[i*nD+j]
        //                << "\T " << creepVel[i*nD+j] / (Units::mm/Units::year)
		<< std::endl;

}

/* end test_loading.cpp */
