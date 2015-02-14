/**
 * Load the fault and compute earthquakes and creep,
 *
 * test of main algorithm including heat interactions
 *
 */
#include "thermomec.h"

using namespace std;

const int nl = 128, nd = 32; // dimensions
const double L = 70*km , D = 17.5*km;

//======================================================================
int main(int argc, char * argv[]) {

  /* Output name of program */
  cout << "Program: " << argv[0] << endl;

  /**
   * Initialization of static properties
   * */
  Halfspace *H = new Halfspace(); // Initialize the halfspace with default values
  Fault F( H , nl , nd , L , D ); // Initialize a fault with default values

  F.initCreepmask_exp(); // use exponential smoothing for horizontal boundary conditions

  vector<double> tmp;
  stressdrop_unif( tmp , nl , nd , 12*bar , 6*bar ); // uniform distribution of stress drops 6-18 bar
  F.initStressDrops(tmp);

  twoLayerActEnergy( tmp , nl , nd , int(0.5+((5*km/D)*nd)) , 100*kJ/mole , 130*kJ/mole ); // two layered activation energy distribution
  F.initActEnergy( tmp );

  /**
   * Set initial values of slip deficit to build some stress on the fault
   * */
  F.initSlipDef("init_slip_def.txt"); // Input file is single column row ordered, units = m,
  //  computed in separate program by inverting the stiffness matrix for desired stress profile
  // Strength profile must have same friction and creep values as initialized above
  //F.init_slipdef( 0.95 , vp ); // TODO: remove need for input file init the slip deficit to 95 % strength

  /**
   * Initializatin of algorithm Constants
   * */
  static const double vp = 35 * mm / year; // plate velocity
  const double slip_velocity = 6 * km /second;     // slip velocity controls time taken for an earthquake


  /**
   * Run algorithm without heat to set up initial distribution of stress
   * */
  double creep_step_for_numerics = 1 * m; // maximum creep distance controls time step
  int maxNhypocenter = 10; // don't allow more than ten hypocenters
  double time0 = 25.0*year;
  double Nyears = 125.0; // Number of years to run the algorithm
  double time_step; // variable time step
  double time = time0;
  int it = 0; // iteration counter
  vector<double> slip_deficit_rate(nl*nd); // rate of loading

  /* Load the fault for time0 years on top of the initial slip deficit*/
  F.load( time0 *vp );

  /* calculate the stress in each cell */
  F.computeStress();

  while ( time <= time0+Nyears*year )
    {

      int hi(0) , hj(0); // hypocenter parameters

      /* If critical define a new earthquake */
      int ncrit = F.is_critical( hi , hj );
      if ( ncrit > 0 ){

	cout << "\nHypocenter at [" << hi << ", " << hj << "], nCrit = " << ncrit << " " ;

	/* Define new earthquake, but don't record in catalog or add to fault for heat*/
	Earthquake this_eq( hi , hj , slip_velocity , time , F );

	cout << "t = " << time/year << " yr"
	     << " / P0 = " << this_eq.getPotency( F )/(cm*km*km) << " km^2 cm\n";

	/* calculate the stress in each cell */
	F.computeStress();

      }

      /* Compute loading rate given creep */
      vector<double> slip_deficit_rate(nl*nd);
      double maxVc = F.getSlipDefRate( vp , slip_deficit_rate );

      //time_step = init_time_step; // constant time step
      time_step =  creep_step_for_numerics / maxVc; // Compute time step to according to allowed slip


      /* check slip deficit rate, gives warning */
      F.checkForUnderOvershoot( vp , time_step , slip_deficit_rate );

      /* Reduce time step so at most maxNhypocenter hypocenters */
      F.adjustTimeStep( maxNhypocenter , slip_deficit_rate , time_step );

      cout << "t = " << time/year << " yr, dt = " << time_step/day << " days\n";

      /* Define Creep rate */
      CreepEvent this_creep( vp , slip_deficit_rate , time_step , time , F );

      /* calculate the stress in each cell*/
      F.computeStress();

      it++;
    }// ---- END algorithm 1

  /* Print the stress and slip deficit on the fault */
  char ofilename0[] = "thermomec.t50.stress-slipdef.out";
  cout << "Writing to " << ofilename0 << "\n";
  ofstream ofile0(ofilename0); // open file

  for( int i = 0 ; i < nl ; i++ )
    for( int j = 0 ; j < nd ; j++ )
      ofile0 << i << " " << j << " "
	     << F.stress(i,j)/MPa <<  " "
	     << F.slipdef(i,j) / m
	     << endl;
  ofile0.close();


  /**
   * Run algorithm with heat and record the data.
   */
  time0 = time;
  creep_step_for_numerics = 1 * cm; // maximum creep distance controls time step
  maxNhypocenter = 1; // don't allow more than one hypocenter
  it = 0; // iteration counter
  Nyears = 10.0; // Number of years to run the algorithm

  /* Catalogs for earthquakes and creep events */
  vector <Earthquake> eqcat;
  vector <CreepEvent> creepcat;

  while ( time <= time0+Nyears*year && it <= 1e6 )
    //while ( it <= 20 )
    {
      /* hypocenter parameters */
      int hi(0) , hj(0);

      /* If critical define a new earthquake */
      int ncrit = F.is_critical( hi , hj );
      if ( ncrit > 0 ){

	cout << "\nHypocenter at [" << hi << ", " << hj << "], nCrit = " << ncrit << " " ;

	/* Define new earthquake */
	Earthquake this_eq( hi , hj , slip_velocity , time , F );

	eqcat.push_back( this_eq ); // add to the catalog

	F.add_eqk_record( this_eq.getCellslips_ptr() , // add record of slip to fault
			  this_eq.getStrikeIdxs_ptr() , this_eq.getDepthIdxs_ptr() );

	cout << "t = " << time/year << " yr / Eq # = " << eqcat.size() << " / P0 = " << this_eq.getPotency( F )/(cm*km*km) << " km^2 cm\n";

	/* calculate the stress in each cell */
	F.computeStress();

	/* calculate the temperature in each cell */
	F.computeTemp( time , 1000.0*year , 1000.0*year );

      }

      /* Compute loading rate given creep */
      double maxVc = F.getSlipDefRate( vp , slip_deficit_rate );

      /* Compute the time step */
      time_step =  creep_step_for_numerics / maxVc; // Compute time step to according to allowed slip

      /* Correct slip deficit rate, gives warning if changed*/
      F.checkForUnderOvershoot( vp , time_step , slip_deficit_rate );

      /* Correct time step so only one hypocenter */
      F.adjustTimeStep( maxNhypocenter , slip_deficit_rate , time_step );

      /* Get the heat rate  by computing the slip deficit rate after half a time step */
      do {
	maxVc = F.recomputeSlipDefRate(vp , time , time_step , 1000*year , 1000*year , slip_deficit_rate  );

	/* check we don't need to reduce the time step for new loading rate*/
	if( creep_step_for_numerics / maxVc  < time_step ) time_step = creep_step_for_numerics / maxVc;

	/* Correct slip deficit rate, gives warning if changed*/
	F.checkForUnderOvershoot( vp , time_step , slip_deficit_rate );

	/* Check for only one hypocenter - redo calculation for smaller time step if so*/
      } while( F.adjustTimeStep( maxNhypocenter , slip_deficit_rate , time_step ) );

      /* Define Creep Event */
      CreepEvent this_creep( vp , slip_deficit_rate , time_step , time , F );

      cout << "t = " << time/year << " yr "
	   << ", dt = " << time_step/day << " day"
	   << "/ Creep P0 = " << this_creep.getPotency( F )/(cm*km*km) << " cm km^2.\n";

      /* Update stress and heat rate */
      F.computeStress();/* calculate the stress in each cell */

      /* Add record of the creep to the fault and to catalog */
      F.add_creep_record( this_creep.getCreepslips_ptr() , time );
      creepcat.push_back( this_creep );

      /* Update temperature */
      F.computeTemp( time , 1000.0*year , 1000.0*year );
      //F.checkTemp( TmQtz );  // ...tmp


      it++;
    } // ---- END algorithm 2

  /* Results to output file */

  /* Print the Eq cat results */
  char ofilename1[] = "thermomec.eqcat.out";
  cout << "Writing to " << ofilename1 << "\n";
  ofstream ofile1( ofilename1 ); // open file

  for( unsigned int i = 0 ; i < eqcat.size() ; i++ ) {
    double x , z; // get hypocenter coordinates
    F.getCoords( eqcat[i].getHi() , eqcat[i].getHj() , x , z );

    ofile1 << eqcat[i].getTime() /year << " " // Time
	   << x / km << " " << z / km << " " // Hypocenter location
	   << eqcat[i].getnFailed()*F.cell_area() / (km*km) << " " // Rupture Area
	   <<  eqcat[i].getPotency( F ) / (cm*km*km) << " " // Scalar Potency
	   << eqcat[i].getMagnitude( F ) << endl; // Magnitude
  }
  ofile1.close();

  /* Print the creep cat results */
  char ofilename2[] = "thermomec.creepcat.out";
  cout << "Writing to " << ofilename2 << "\n";
  ofstream ofile2( ofilename2 ); // open file

  for( unsigned int i = 0 ; i < creepcat.size() ; i++ ) {

    ofile2 << creepcat[i].getTime()/year << " " // Time
	   << creepcat[i].getnFailed()*F.cell_area() / (km*km) << " " // Rupture Area
	   <<  creepcat[i].getPotency( F ) / (cm*km*km) << " " // Scalar Potency
	   << creepcat[i].getMagnitude( F ) << " "
	   << creepcat[i].getAvgHeatRate( F.cell_area() ) / (J/second) << endl; // Magnitude
  }
  ofile2.close();

  /* Print the temperature results */
  char ofilename3[] = "thermomec.temp.out";
  cout << "Writing to " << ofilename3 << "\n";
  ofstream ofile3(ofilename3); // open file

  for( int i = 0 ; i < nl ; i++ )
    for( int j = 0 ; j < nd ; j++ )
      ofile3 << i << " " << j << " "
	     << ( F.taus(i,j) - F.taua(i,j) )/MPa <<  " "
	     << F.temperature(i,j) - 273 << " "
	     << F.temperature(i,j) - F.background_temperature(i,j)
	     << endl;
  ofile3.close();

  /* Print the slip partition results */
  char ofilename4[] = "thermomec.partition.out";
  cout << "Writing to " << ofilename4 << "\n";
  ofstream ofile4(ofilename4); // open file

  for( int i = 0 ; i < nl ; i++ )
    for( int j = 0 ; j < nd ; j++ )
      ofile4 << i << " " << j << " "
	     << F.getStrength(i,j,1.0,vp)/MPa <<  " "
	     << F.getTotalSeisSlip(i,j)/m << " "
	     << F.getTotalCreepSlip(i,j)
	     << endl;
  ofile4.close();
}
//======================================================================
//######################################################################
//######################################################################
