/**
 * test_circlerupt.cpp ---
 *
 * Generate a circular earthquake and output the stress, temp and slip rate
 *
 * Filename: test_circlerupt.cpp
 * Author: Iain W. Bailey
 * Created: Sun Jul 13 21:36:51 2014 (-0400)
 * Version: 1
 * Last-Updated: Sun Jul 20 21:27:20 2014 (-0400)
 *           By: Iain W. Bailey
 *     Update #: 162
 *

 * Change Log:

*/


#include <vector>
#include <iostream>
#include <stdio.h>

#include "units.h"
#include "fileio.h"
#include "temperature_profile.h"
#include "strength_profile.h"
#include "side_creep.h"
#include "creep_law.h"
#include "stiffness_matrix.h"
#include "thermofault.h"
#include "earthquake.h"

#include "parameters.h"

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

  /* Background temperature */
  BackgroundTemperatureProfile T( In::Tsurface, In::Tgradient );

  /* Get the  background temperature */
  std::vector<double> faultBkgdT = T.faultgrid_co( In::nL, In::nD, In::cellHeight );

  /* Compute the stiffness matrix */
  StiffnessMatrix K( In::cellLength, In::nL, In::cellHeight, In::nD, In::rigidity );

  /* Static strength profile */
  StaticStrengthProfile StatStrength( In::tau0, In::fs, In::dsigmadz );

  /* Set up the creep parameters over entire fault including creep mask*/
  int nCells = In::nL*In::nD;
  std::vector<double> faultE( nCells, In::activationEnergy);
  std::vector<double> faultA( nCells, In::arrhAmplitude);
  std::vector<double> faultn( nCells, In::stressExponent );
  sidecreepmask( faultA, faultE, faultn, faultBkgdT,
                 In::nL, In::nD, In::cellLength, In::xDB );

  /* Set up the strength drops */
  std::vector<double> strengthDrops( nCells, In::dynStrengthDrop );

  /*Initialize fault array properties */
  std::vector<double>
    bkgdT(nCells),
    initStress(nCells),
    taus(nCells),
    taud(nCells);

  /* Loop depth */
  for( int j=0; j<In::nD; j++ )
    {
      double depth = ((double)j + 0.5)*In::cellHeight;

      /* Loop along strike */
      for( int i=0; i<In::nL; i++ )
	{
          double xPos = ((double)i + 0.5)*In::cellHeight;

          /* Cell index (column ordered) */
          int iCell = i*In::nD+j;

          /* Set the static strength */
	  taus[iCell] = StatStrength(depth);

          /* Dynamic strength */
	  taud[iCell] = taus[iCell] - strengthDrops[iCell];

          /* Creep law at this cell */
	  CreepLaw C( faultA[iCell], faultE[iCell], faultn[iCell] );

          /* Init stress so earthquake is critical */
          if( pow( xPos - In::hx, 2 ) + pow( depth - In::hz, 2 ) <=
              pow( In::radiusRupture, 2 ) )
            /* Set above static strength within circular rupture area */
            initStress[iCell] =  1.1*taus[iCell];
          else
            /* Set below static or dynamic strength outside of rupture */
            initStress[iCell] = std::min( 0.9*StatStrength(depth),
                                          C.stress( In::loadingStrainRate, faultBkgdT[iCell] ));
	} //--END strike loop
    } //-END depth loop

  /* Set up the cooling model */
  FaultCoolingModel *FCM = new FaultCoolingModel;

  /** Initialize the fault */
  ThermoFault F( In::nL, In::nD, In::cellLength, In::cellHeight, In::faultWidth,
                 taus, taud, In::dosCoef, faultA, faultn, faultE, K,
                 initStress, faultBkgdT, FCM );

  /* Set containers for storing the fault info */
  std::vector<double>
    stress(nCells),
    temperature(nCells),
    creepVel(nCells);

  /* Containers for recording the slip */
  std::vector<double> totalCreepSlip( nCells, 0.0), totalEqkSlip( nCells, 0.0);

  /* Save the initial stress */
  faultvec2txtfile( In::ofilename_initStress, initStress*(1.0/Units::MPa), In::nL, In::nD );
  faultvec2txtfile( In::ofilename_initTemp, faultBkgdT, In::nL, In::nD );

  /* Calculate the stress and temperature*/
  double time = 0.0;
  F.getStress( stress );
  F.getTemperature( temperature, time );

  /* Print header for the earthquake catalog */
  std::cout << "Time_yr, x_km, z_km, Mag_P, Mag_W, Area_km2, StressDrop_MPa\n";

   /** Run the algorithm */
  int i = 0;
  while( i < In::nTimeMax && time < In::maxTime )
    {
      /* Get the creep velocity on the fault */
      F.getCreepVelocity( creepVel, stress, temperature);

      /* Compute the time to failure */
      double timeStep = F.estimateTimeToFailure( stress, creepVel, In::plateVelocity );

      /* Check and adjust the time step */
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

      /* Get the new stress and temperature*/
      F.getStress( stress );

      /* Calculate whether there are any hypocenters */
      unsigned int iHypo, jHypo;
      unsigned int nCrit = F.nCriticalCells(iHypo, jHypo, stress );

      /* Warning message for multiple hypocenters, otherwise do nothing */
      if( nCrit > 1 ) std::cerr << "\nWARNING: Multiple hypocenters (" << nCrit << ") ";

      /* Compute earthquake if there are hypocenters */
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
        std::vector<double> eqkHeat( nCells, 0.0 );
        F.computeEarthquake( stress, time, iHypo, jHypo, In::slipVelocity,
                             eqkSlip, eqkHeat );

        /* Compute the new slip deficit */
        std::vector<double> slipDeficitPost;
        F.getSlipDeficit( slipDeficitPost );

        /* Compute the earthquake properties */
        Earthquake thisEQ( iHypo, jHypo, timePre, In::cellLength*In::cellHeight,
                       slipDeficitPre, slipDeficitPost, stressPre, stress );

        /* Record the total earthquake slip */
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
      } /*END earthquake */

      /* Compute the new temperature */
      F.getTemperature( temperature, time );

      /* Update iteration count */
      i++;

      /* Add new line to stderr */
      std::cerr << std::endl;

    } /*END fault algorithm */

  /* Break the line */
  std::cerr << std::endl;

  /* Report why we finished the loop */
  if( i == In::nTimeMax ) std::cerr << "Max number of iterations reached\n";
  else std::cerr << "Max time reached\n";

  /* Output files */
  faultvec2txtfile( In::ofilename_creepSlip, totalCreepSlip, In::nL, In::nD ); // total creep slip
  faultvec2txtfile( In::ofilename_eqSlip, totalEqkSlip, In::nL, In::nD ); // total eq slip

  std::vector<double> slipDeficit;
  F.getSlipDeficit( slipDeficit );
  faultvec2txtfile( In::ofilename_slipdef, slipDeficit*(1.0/Units::m), In::nL, In::nD );

  /* Save the final stress */
  faultvec2txtfile( In::ofilename_finalStress, stress*(1.0/Units::MPa), In::nL, In::nD );
  faultvec2txtfile( In::ofilename_finalTemp, temperature, In::nL, In::nD );

 // /* Init counters */
  // double time=0.0;
  // int it = 0;

  // /* Initialize the fault */
  // Fault F( H , nl, nd, L, D, thickness, xDB, cohesion, fs, doc, vp,
  //          stressdrops, actEnergies);

  // /* Set background stress */
  // F.initBkdStress( 0.99, vp );
  // F.computeStress();

  // /* Compute temp */
  // F.computeTemp( time );

  // write_stress( it, string(argv[0]), time, F ); // write stress
  // write_sliprate( it, string(argv[0]), time, F ); // write slip rate
  // write_temp( it, string(argv[0]), time, F ); // write temp

  // /* Make circular patch above static strength*/
  // vector<double> stress0;
  // double x, y; // fault positions
  // double r2 = r*r; // radius squared

  // for( int i = 0 ; i < nl ; i ++ )
  //   for ( int j = 0 ; j < nd ; j ++ ){

  //     F.getCoords( i, j, x, y );
  //     if( ( x - hx )*( x - hx ) + ( y - hy )*( y - hy ) <= r2 )
  //       stress0.push_back( 1.01*F.taus(i,j) );
  //     else
  //       stress0.push_back( F.stress(i,j) );

  //     }

  // F.initBkdStress( stress0 );   /* Change the background stress */
  // F.computeStress();
  // it++;
  // write_stress( it, string(argv[0]), time, F ); // write stress
  // write_temp( it, string(argv[0]), time, F ); // write temp
  // write_sliprate( it, string(argv[0]), time, F ); // write slip rate

  // /* Compute earthquake */
  // int hi(0) , hj(0), ncrit(0);     // hypocenter parameters
  // vector<Earthquake> eqcat;

  // ncrit = F.is_critical( hi , hj );
  // if ( ncrit > 0 ){

  //   cout << "\nHypocenter at [" << hi << ", " << hj
  //        << "], nCrit = " << ncrit << "\n" ;

  //   /* Define new earthquake */
  //   Earthquake this_eq( hi , hj , time , F );
  //   eqcat.push_back( this_eq );

  //   /*update time */
  //   time += this_eq.getMaxSlip()/slip_velocity;

  //   F.add_eqk_record( this_eq.getCellslips_ptr() , // add record of slip to fault
  //       	      this_eq.getStrikeIdxs_ptr() ,
  //       	      this_eq.getDepthIdxs_ptr() );

  //   double xhypo, yhypo;
  //   F.getCoords( hi, hj, xhypo, yhypo );
  //   cout << "t = " << time/year << " yr"
  //        << " / Hypo = [ " << xhypo/km << "," << yhypo/km << " ] km"
  //        << " / A = " << this_eq.getnFailed()*F.cell_area() / (km*km) << " km^2"
  //        << " / M = " << this_eq.getMagnitude( F )
  //        << " / P0 = " << this_eq.getPotency( F )/(cm*km*km) << " km^2 cm\n";

  // }// --- END EQK

  // /* Compute temp */
  // F.computeTemp( time );

  // /* Write files */
  // it++;
  // write_stress( it, string(argv[0]), time, F );
  // write_temp( it, string(argv[0]), time, F ); // write temp
  // write_sliprate( it, string(argv[0]), time, F ); // write slip rate

  // /* calculate temperature in each cell */
  // F.computeTemp( time );

  // //  Earthquake this_eqk = impose_earthquake( F , time );
  // //F.computeStress();


}


//
// test_circlerupt.cpp ends here
