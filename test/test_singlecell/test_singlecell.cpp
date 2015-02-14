/**
 * Make a fault with just one cell.
 *
 * Impose an earthquake then compute temperature and creep in the cell
 * afterwards.
 *
 *
 */
#include <vector>
#include <iostream>
#include <stdexcept> // for runtime_error

#include "input_parameters.h" // Parameters used to set up the fault

#include "temperature_profile.h"
#include "strength_profile.h"
#include "side_creep.h"
#include "creep_law.h"
#include "stiffness_matrix.h"
#include "thermofault.h"
#include "fault_cooling.h"

//------------------------------------------------------------------------------
int main(int argc, char * argv[]) {

  std::cerr << "Program: " << argv[0] << " No. Args: " << argc << std::endl;

  /* Stiffness matrix */
  StiffnessMatrix K( In::cellLength, 1, In::cellDepth, In::cellHeight, 1,
		     30*Units::GPa ); // rigidity of the halfspace

  /* Background temperature */
  BackgroundTemperatureProfile T( In::Tsurface, In::Tgradient );

  /* Static strength profile */
  StaticStrengthProfile StatStrength( In::tau0, In::fs, In::dsigmadz );

  /* Static strength of the fault */
  std::vector<double> taus( 1, StatStrength( In::cellDepth ) );
  std::vector<double> taud( 1, taus[0] - In::dynStrengthDrop );
  std::vector<double> bkgdT( 1, T(In::cellDepth) );

  /* Creep parameters of fault */
  std::vector<double> faultA( 1, In::arrhAmplitude );
  std::vector<double> faultE( 1, In::activationEnergy );
  std::vector<double> faultn( 1, In::stressExponent );

  /** Initialize the fault */
  ThermoFault F( 1, 1, In::cellLength, In::cellHeight, In::faultWidth,
                 taus, taud, In::dosCoef,
                 faultA, faultn, faultE,
                 K,
                 taus,
                 bkgdT,
                 new FaultCoolingModel() );

  /* Containers for the fault stress and temperature */
  std::vector<double> stress(1, 0.0);
  std::vector<double> temperature(1, 0.0);
  std::vector<double> creepVel(1,0.0);

  /* Containers for the hypocenter */
  unsigned int iHypo, jHypo;

  /* Output headers*/
  std::cout << "time_yr, stress_MPa, temperature_K, creepVelocity_mm/yr\n";

  /* Set the time */
  double time = 0.0;
  while( time < In::maxTime )
    {
      /* Get the stress and temperature on the fault */
      F.getStress( stress );
      F.getTemperature( temperature, time );

      /* get the creep velocity */
      F.getCreepVelocity( creepVel, stress, temperature );

     /* Output */
      std::cout << time/Units::year << " "
                << stress[0]/Units::MPa << " "
                << temperature[0] << " "
                << creepVel[0]/(Units::mm/Units::year) << std::endl;

      /* Check if any earthquakes */
      int nHypo = F.nCriticalCells( iHypo, jHypo, stress );
      if( nHypo > 0 )
	{
          std::cerr << "Earthquake at time " << time/Units::day << " days\n";

	  /* Compute earthquake */
          std::vector<double> eqkSlip(1, 0.0), eqkHeat(1, 0.0);
          F.computeEarthquake( stress, time, iHypo, jHypo, In::slipVelocity,
                               eqkSlip, eqkHeat );

	  /* Get temperature on the fault */
	  F.getTemperature( temperature, time );
          F.getCreepVelocity( creepVel, stress, temperature );

          /* Output */
          std::cout << time/Units::year << ", "
                    << stress[0]/Units::MPa << ", "
                    << temperature[0] << ", "
                    << creepVel[0]/(Units::mm/Units::year) << std::endl;
	}

      /* Compute loading time, minimum time to failure or time for max slip*/
      double timeToFail = F.estimateTimeToFailure( stress, creepVel,
                                                   In::plateVelocity );

      double maxVel = *(std::max_element( creepVel.begin(), creepVel.end() ));
      double timeForSlip = In::maxCreepSlip / maxVel;
      double timeStep = std::min( std::min( timeToFail, timeForSlip ),
                                  In::maxTimeStep );

      // std::cout << "\tCreep velocity: " << maxVel/(Units::mm/Units::year)
      //           << " mm/yr\n"
      //           << "\ttime to fail: " << timeToFail/Units::day << " days, "
      //           << "time for slip: " << timeForSlip/Units::day << " days\n"
      //           << "\ttimeStep = " << timeStep/Units::day << " days.\n";

      if( timeStep <=0 )  throw std::runtime_error("Non-positive stress");

      /* Apply loading */
      F.loadFault( In::plateVelocity, creepVel, timeStep );

      /* Update time */
      time += timeStep;

    }

  // //const double tInit = 15*year; // initial loading time
  // const double slip_velocity = 6 * km /second;     // slip velocity
  // const int maxNhypocenter = 1; // don't allow more than one hypocenter
  // const double creep_step_for_numerics = 35 * mm; // creep distance constant

  // vector <Earthquake> eqcat; // catalog for recording the eqks
  // vector <CreepEvent> creepcat; // catalog for recording the creep

  // double x , y ; // coordinates
  // F.getCoords(0,0,x,y);
  // cout << "Single Cell created.  Central Coordinates: x = " << x/km << " km , y = " << y/km << " km\n";
  // cout << "Cell Area = " << F.cell_area() /(km*km) << " km^2\n";
  // cout << "taus = " << F.taus(0,0)/MPa << " MPa, "
  //      << "taud = " << F.taud(0,0)/MPa << " MPa, "
  //      << "taua = " << F.taud(0,0)/MPa << " MPa, "
  //      << "taub = " << F.taub(0,0)/MPa << " MPa\n";

  // /* Compute stress */
  // F.computeStress();
  // cout << "tau = " << F.stress(0,0)/MPa << " MPa\n";

  // /* predict time to failure */
  // double tInit = F.predictTimeToFailure(0,0,vp);
  // double time = tInit;
  // double eq_cutoff_time = 1000*year, creep_cutoff_time =1000*year;

  // F.load( vp*tInit ); // load fault

  // /* Compute stress and temp */
  // F.computeStress();
  // F.computeTemp( tInit , eq_cutoff_time , creep_cutoff_time );

  // cout << "Loaded with u = " << vp*tInit/m << " m over " << tInit/year << " yrs\n";
  // cout << "tau = " << F.stress(0,0)/MPa << " MPa\n";

  // double time_step; // time step ... variable
  // vector<double> slip_deficit_rate(nl*nd); // rate of loading
  // vector<double> heat_rate( nl*nd ); // rate of heating
  // int it = 0; // iteration counter

  // //while( time <= tInit + nYears*year )
  // while( it < 3 )
  //   {
  //     cout << "\nIteration: " << it << endl;
  //     cout << "t = " << time/year << " yr / "
  // 	   << "udef = " << F.slipdef(0,0)/m << " m / "
  // 	   << "tau = " << F.stress(0,0)/MPa << " MPa / "
  // 	   << " T = " <<  F.temperature(0,0) - 273 << " ^oC\n";

  //     // check for earthquake
  //     int hi(0) , hj(0);

  //     /* If critical define a new earthquake */
  //     int ncrit = F.is_critical( hi , hj );
  //     if ( ncrit > 0 ){

  // 	cout << "\nHypocenter at [" << hi << ", " << hj << "], nCrit = " << ncrit << " " ;

  // 	/* Define new earthquake */
  // 	Earthquake this_eq( hi , hj , time , F );

  // 	/* Update time */
  // 	time += this_eq.getMaxSlip() / slip_velocity;

  // 	eqcat.push_back( this_eq ); // add to the catalog

  // 	F.add_eqk_record( this_eq.getCellslips_ptr() , // add record of slip to fault
  // 			  this_eq.getStrikeIdxs_ptr() , this_eq.getDepthIdxs_ptr() );

  // 	cout << "Eq # = " << eqcat.size()
  // 	     << ", ML = " << this_eq.getMagnitude( F ) << ", "
  // 	     << " du = " << this_eq.getMaxSlip()/m << " m\n";

  // 	/* calculate the stress in each cell */
  // 	F.computeStress();

  // 	/* calculate the temperature in each cell */
  // 	F.computeTemp( time , eq_cutoff_time , creep_cutoff_time );

  // 	cout << "t = " << time/year << " yr / "
  // 	     << "udef = " << F.slipdef(0,0)/m << " m / "
  // 	     << "tau = " << F.stress(0,0)/MPa << " MPa / "
  // 	     << " T = " <<  F.temperature(0,0) - 273 << " ^oC\n\n";
  //     }

  //     /* Compute loading rate given creep */
  //     double maxVc = F.getSlipDefRate( vp , slip_deficit_rate );
  //     cout << "Slip rate based on current temp and stress = " << maxVc /(mm/year) << " mm/yr\n";

  //     /* Compute the time step */
  //     time_step =  creep_step_for_numerics / maxVc; // Compute time step to according to allowed slip
  //     cout << "Time step for max slip allowed = " << time_step /year << " year\n";

  //     time_step = min(time_step , F.predictTimeToFailure(0,0,vp-maxVc) );

  //     /* Correct time step so only one hypocenter */
  //     //F.adjustTimeStep( maxNhypocenter , slip_deficit_rate , time_step );

  //     cout << "Time step after check against failure stress = " << time_step /year << " year\n";

  //     /* Get the heat rate  by computing the slip deficit rate after half a time step */
  //     do {
  // 	maxVc = F.recomputeSlipDefRate(vp , time , time_step ,
  // 				       eq_cutoff_time , creep_cutoff_time ,
  // 				       slip_deficit_rate );

  // 	/* check we don't need to reduce the time step for new loading rate*/
  // 	if( creep_step_for_numerics / maxVc  < time_step ) time_step = creep_step_for_numerics / maxVc;

  // 	/* Check for only one hypocenter - redo calculation for smaller time step if so*/
  //     } while( F.adjustTimeStep( maxNhypocenter , slip_deficit_rate , time_step ) );
  //     cout << "Time step accounting for heat generation while creeping = "
  // 	   << time_step /year << " year\n";
  //     cout << "Slip rate accounting for heat generation while creeping = "
  // 	   << maxVc /(mm/year) << " mm/yr\n";

  //     /* Define Creep Event */
  //     CreepEvent this_creep( vp , slip_deficit_rate , time_step , time , F );

  //     /* Update time */
  //     time += time_step;

  //     cout << "t = " << time/year << " yr "
  // 	   << ", Creep ML = " << this_creep.getMagnitude( F ) << ", "
  // 	   << "du = " << this_creep.getPotency( F )/F.cell_area()/m << " m.\n";

  //     /* Update stress and heat rate */
  //     F.computeStress();/* calculate the stress in each cell */
  //     F.getHeatRate( vp , slip_deficit_rate , heat_rate );/* calculate the heat rate at end of time step */

  //     /* Add record of the creep to the fault and to catalog */
  //     F.add_creep_record( this_creep.getCreepslips_ptr()  , time );
  //     creepcat.push_back( this_creep );

  //     /* Update temperature */
  //     F.computeTemp( time , eq_cutoff_time , creep_cutoff_time );

  //     it++;
  //   } // ---- END algorithm


}
//======================================================================
//######################################################################
//######################################################################
