/**
 * test_bz1996.cpp ---
 *
 *  Test a fault set up that is analogous to the BZ 1996 model,
 *  i.e. stress-dependent, but no temperature dependent creep.  Note it is not
 *  exactly the same because the creep is set up in a different way
 *
 *
 * Filename: test_bz1996.cpp
 * Description:
 * Author: Iain W. Bailey
 * Created: Sun Apr  6 14:15:22 2014 (-0400)
 * Version: 1
 * Last-Updated: Sun Aug 24 15:46:13 2014 (-0400)
 *           By: Iain W. Bailey
 *     Update #: 274
 *

 * Change Log:
 *
 *
 *
 */
#include <vector>
#include <iostream>
#include <string.h>
#include <stdio.h>

/* Parameters used in this program */
#include "input_parameters.h"

/* Functions for the calculations */
#include "fileio.h"
#include "temperature_profile.h"
#include "strength_profile.h"
#include "benzion_creep.h"
#include "creep_law.h"
#include "stiffness_matrix.h"
#include "fault.h"
#include "earthquake.h"

/** Make operator for multiplying a vector by constant scalar */
std::vector<double> operator*( const std::vector<double>& v, double alpha)
{
  std::vector<double> v2 = v;
  for( unsigned int i = 0; i<v2.size(); i++ ) v2[i] *= alpha;
  return v2;
}

//------------------------------------------------------------------------------
int main(int argc, char * argv[])
{

  std::cerr << "Program: " << argv[0] << " No. Args: " << argc << std::endl;

  /* Background temperature */
  BackgroundTemperatureProfile T( In::Tsurface, In::Tgradient);

  /* Static strength profile */
  StaticStrengthProfile StatStrength( In::tau0, In::fs, In::dsigmadz );

  /* Dynamic strength drop */
  int nCells = In::nL*In::nD;
  std::vector<double> strengthDrops( nCells );
  if( !txtfile2faultvec( strengthDrops, In::ifile_strengthdrops, In::nL, In::nD ) ){
    return -1;
  }

  /* Compute the stiffness matrix */
  StiffnessMatrix K( In::cellLength, In::nL, In::cellHeight, In::nD, In::rigidity );

  /* Calculate the A value according to the  creep mask used in BZ1996 */
  std::vector<double> faultA(nCells);
  faultcreep_bz96( faultA, In::nL, In::nD, In::cellLength, In::cellHeight,
                   In::xBD, In::xBD,
                   In::loadingStrainRate, StatStrength( In::xBD ) );

  /* Set up the creep parameters over entire fault including creep mask*/
  std::vector<double> faultE(nCells, 0.0);
  std::vector<double> faultn(nCells, 3.0 );

  /*Initialize fault array properties */
  std::vector<double> bkgdT(nCells), initStress(nCells), taus(nCells), taud(nCells);
  for( int j=0; j<In::nD; j++ )
    {
      double depth = ((double)j + 0.5)*In::cellHeight;
      for( int i=0; i<In::nL; i++ )
	{
          int iCell = i*In::nD+j;
	  taus[iCell] = StatStrength(depth);
	  taud[iCell] = taus[iCell] - strengthDrops[iCell];

          /* Calculate the creep strength and the strength envelope */
          CreepLaw C0( faultA[iCell], 0.0, 3.0 );
	  initStress[iCell] = std::min( StatStrength(depth),
                                        C0.stress( In::loadingStrainRate, bkgdT[iCell] ));
	}
    }

  /* Save the initial stress */
  const std::string ofilename_initStress( argv[0] + In::ofilesuffix_initStress );
  faultvec2txtfile( ofilename_initStress, initStress*(1.0/Units::MPa), In::nL, In::nD );

  /** Initialize the fault */
  Fault F( In::nL, In::nD, In::cellLength, In::cellHeight, In::faultWidth,
           taus, taud, In::dosCoef, faultA, faultn, faultE, K, initStress );

  /* Set containers for storing the fault info */
  std::vector<double> stress(nCells), creepVel(nCells);

  /* Containers for recording the slip */
  std::vector<double> totalCreepSlip( nCells, 0.0), totalEqkSlip( nCells, 0.0);

  /* Calculate the stress */
  F.getStress( stress );

  /* Print header for the earthquake catalog */
  std::cout << "Time_yr, x_km, z_km, Mag_P, Mag_W, Area_km2, StressDrop_MPa\n";

  /** Run the algorithm */
  int i = 0;
  double time = 0.0;
  while( i<In::nTimeMax && time < In::maxTime )
    {
      /* Get the creep velocity on the fault */
      F.getCreepVelocity( creepVel, stress, bkgdT);

      /* Compute the time to failure */
      double timeStep = F.estimateTimeToFailure( stress, creepVel, In::plateVelocity );

      /* Adjust the time step */
      if( timeStep < 0 ) {
        /* Negative implies something went wrong */
        std::cerr << "ERROR: negative time step\n";
        return -1;
      } else if( timeStep > In::maxTimeStep ){
        /* Don't let it get too big or creep rates will be inaccurate */
        std::cerr << "Time " << time/Units::year << ": Using maximum time step. ";
        timeStep = In::maxTimeStep;
      } else if( timeStep < In::minTimeStep ){
        /* Don't let it get too small or we will be waiting forever */
        std::cerr << "Time " << time/Units::year << ": Using minimum time step. ";
        timeStep = In::minTimeStep;
      } else{
        std::cerr << "Time " << time/Units::year;
      }

      /* Load the fault */
      F.loadFault( In::plateVelocity, creepVel, timeStep );

      /* Add on the creep slip */
      for( unsigned int k=0; k<nCells; k++ )
        totalCreepSlip[k] += creepVel[k]*timeStep;

      /* Update the time */
      time += timeStep;

      /* Get the new stress */
      F.getStress( stress );

      /* Calculate whether there are any hypocenters */
      unsigned int iHypo, jHypo;
      unsigned int nCrit = F.nCriticalCells(iHypo, jHypo, stress );
      if( nCrit > 1 ) std::cerr << "\nWARNING: Multiple hypocenters (" << nCrit << ") ";
      if( nCrit > 0 ){
        /* Output a dot for progress to terminal */
        std::cerr << " *EQ* ";

        /* Compute the slip deficit before the earthquake */
        std::vector<double> slipDeficitPre;
        F.getSlipDeficit( slipDeficitPre );

        /* Store the prior stress */
        std::vector<double> stressPre = stress;

        /* Store the prior time */
        double timePre = time;

        /* Compute the earthquake */
        std::vector<double> eqkSlip( nCells, 0.0 );
        F.computeEarthquake( stress, time, iHypo, jHypo, In::slipVelocity, eqkSlip );

        /* Compute the new slip deficit */
        std::vector<double> slipDeficitPost;
        F.getSlipDeficit( slipDeficitPost );

        /* Compute the earthquake properties */
        Earthquake thisEQ( iHypo, jHypo, timePre, In::cellLength*In::cellHeight,
                       slipDeficitPre, slipDeficitPost, stressPre, stress );

        /* Record the total earthquake slip */
        if( i > 0 ){
          for( unsigned int k=0; k<nCells; k++ )
            totalEqkSlip[k] += slipDeficitPre[k] - slipDeficitPost[k];

          /* Output the earthquake catalog: time, x, z, magP, magW, area, stressdrop*/
          printf( "%9.6E, %6.2f, %5.2f, %4.2f, %4.2f, %6.2f, %5.2f\n",
                  timePre/Units::year,
                  (iHypo+0.5)*In::cellLength/Units::km ,
                  (jHypo+0.5)*In::cellHeight/Units::km ,
                  thisEQ.potMagnitude() ,
                  thisEQ.momMagnitude( In::rigidity ) ,
                  thisEQ.ruptureArea() / (Units::km*Units::km) ,
                  thisEQ.staticStressDrop()/Units::MPa );
        }
      } /*END earthquake */

      i++;
      std::cerr << std::endl;

    } /*END fault algorithm */

  /* Break the line */
  std::cerr << std::endl;

  /* Report why we finished the loop */
  if( i == In::nTimeMax ) std::cerr << "Max number of iterations reached\n";
  else std::cerr << "Max time reached\n";

  /* Output files */
  const std::string ofilename_creepSlip( argv[0] + In::ofilesuffix_creepSlip );
  faultvec2txtfile( ofilename_creepSlip, totalCreepSlip, In::nL, In::nD ); // total creep slip

  const std::string ofilename_eqSlip( argv[0] + In::ofilesuffix_eqSlip );
  faultvec2txtfile( ofilename_eqSlip, totalEqkSlip, In::nL, In::nD ); // total eq slip

  std::vector<double> slipDeficit;
  F.getSlipDeficit( slipDeficit );
  const std::string ofilename_slipdef( argv[0] + In::ofilesuffix_slipdef );
  faultvec2txtfile( ofilename_slipdef, slipDeficit*(1.0/Units::m), In::nL, In::nD );

  /* Save the final stress */
  const std::string ofilename_finalStress( argv[0] + In::ofilesuffix_finalStress );
  faultvec2txtfile( ofilename_finalStress, stress*(1.0/Units::MPa), In::nL, In::nD );

}

/* END */
