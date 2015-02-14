/**
 * Impose an earthquake on the fault and look at the change in heat at
 * various times afterwards
 *
 *
 */
#include "thermomec.h"

using namespace std;

/* Constants for halfspace */
const double rigidity( 30 * GPa ),// rigidity of half space
  amp_arrh( 6.31e-20 *1/(Pa*Pa*Pa*second) ), 
  density( 2.6*g/(cm*cm*cm) ), // PREM density of crustal rock kg/km^3
  specHeat( 790*J/(kg*K) ), // specific heat of granite J/(kg*K)
  diffusivity( 1.E-2*cm*cm/second ),  // diffusivity of rock 10-2 cm^2/s
  coolDist( 100*m ),  // distance where normal temp maintained
  Tsurface( ( 273.15 + 0 ) * K ), // temperature of 283 K at surface
  Tgrad( 20 * K/km ), // degrees K per km gradient
  sigmaNgrad( 280*(bar/km) ) , // rate of normal stress increase
  porePgrad( 100*(bar/km) ) ,
  AE0( 130 *kJ/mole ); // rate of pore pressure increase

/* Constants for fault */
const int nl = 128, nd = 32; // number of grid cells
const double L = 70*km, D = 17.5*km, thickness( 10*cm ), // Fault dimensions
  xDB( 10.0*km ) ,// creeping boundary width
  cohesion( 10*MPa ) , // cohesion & max stress drop
  fs( 0.75 ) , // static coeff of friction
  doc( 1.25 );// dynamic overshoot coeff

const double vp( 35 * mm / year ) // plate velocity 
  , slip_velocity( 4*km/second ) // slip velocity 
  , dtmax( 1*day ) // max time step
  , maxcreepslip( 1*cm );

/* Earthquake parameters */
const double r( 3*km ), hx( 35.5*km ), hy( 5*km ),
  dtau( 5*MPa );

string idxtime_fname("./txtfiles/idxtime_map.txt");

//######################################################################
int main(int argc, char * argv[]) {

  cerr << "Program: " << argv[0] << " No. Args: " << argc << endl;

  //Halfspace *H = new Halfspace(); // Initialize the halfspace 
  Halfspace *H = new Halfspace( rigidity, amp_arrh, density, specHeat, 
				diffusivity, coolDist, Tsurface, Tgrad, 
				sigmaNgrad, porePgrad ); 


  vector<double> stressdrops( nl*nd , dtau ); // taus - taua
  //stressdrop_unif( stressdrops , nl , nd , "unif_0to1.txt", 20*bar , 10*bar ); 
  vector<double> actEnergies( nl*nd, AE0 ); // constant act energy

  /* Init counters */
  double time=0.0;
  int it = 0;

  /* Initialize the fault */
  Fault F( H , nl, nd, L, D, thickness, xDB, cohesion, fs, doc, vp, 
	   stressdrops, actEnergies); 

  /* Set background stress */
  F.initBkdStress( 0.99, vp );
  F.computeStress();

  /* Compute temp */
  F.computeTemp( time );

  /* Write files */
  ofstream ofile;
  ofile.open( idxtime_fname.c_str() );
  write_stress( it, string(argv[0]), time, F ); // write stress
  write_sliprate( it, string(argv[0]), time, F ); // write slip rate
  write_temp( it, string(argv[0]), time, F ); // write temp
  ofile << it << " " << time/year << endl;

  /* Make circular patch above static strength*/
  vector<double> stress0;
  double x, y; // fault positions
  double r2 = r*r; // radius squared

  for( int i = 0 ; i < nl ; i ++ )
    for ( int j = 0 ; j < nd ; j ++ ){
      
      F.getCoords( i, j, x, y );
      if( ( x - hx )*( x - hx ) + ( y - hy )*( y - hy ) <= r2 )
	stress0.push_back( 1.01*F.taus(i,j) );
      else
	stress0.push_back( F.stress(i,j) );

      }

  F.initBkdStress( stress0 );   /* Change the background stress */
  F.computeStress();
  it++;
  write_stress( it, string(argv[0]), time, F ); // write stress
  write_temp( it, string(argv[0]), time, F ); // write temp
  write_sliprate( it, string(argv[0]), time, F ); // write slip rate
  ofile << it << " " << time/hour << endl;

  /* Compute earthquake */
  int hi(0) , hj(0), ncrit(0);     // hypocenter parameters 
  vector<Earthquake> eqcat;

  ncrit = F.is_critical( hi , hj );
  if ( ncrit > 0 ){
    
    cout << "\nHypocenter at [" << hi << ", " << hj 
	 << "], nCrit = " << ncrit << "\n" ;
    
    /* Define new earthquake */
    Earthquake this_eq( hi , hj , time , F );
    eqcat.push_back( this_eq );

    /*update time */
    time += this_eq.getMaxSlip()/slip_velocity;
    
    F.add_eqk_record( this_eq.getCellslips_ptr() , // add record of slip to fault
		      this_eq.getStrikeIdxs_ptr() , 
		      this_eq.getDepthIdxs_ptr() );

    double xhypo, yhypo;
    F.getCoords( hi, hj, xhypo, yhypo );
    cout << "t = " << time/year << " yr"
	 << " / Hypo = [ " << xhypo/km << "," << yhypo/km << " ] km"
	 << " / A = " << this_eq.getnFailed()*F.cell_area() / (km*km) << " km^2"
	 << " / M = " << this_eq.getMagnitude( F ) 
	 << " / P0 = " << this_eq.getPotency( F )/(cm*km*km) << " km^2 cm\n"; 

  }// --- END EQK

  /* Compute temp */
  F.computeTemp( time );

  /* Write files */
  it++;
  write_stress( it, string(argv[0]), time, F ); 
  write_temp( it, string(argv[0]), time, F ); // write temp

  /* Compute slip rate */
  double time_step = dtmax;
  vector<double> slip_deficit_rate;
  time_step = F.getSlipDefRate( slip_deficit_rate, vp, 
				time_step, maxcreepslip, 
				time , true );
  write_sliprate( it, string(argv[0]), time, F , slip_deficit_rate , vp ); // write slip rate
  ofile << it << " " << time/hour << endl;
  cout << "dt = " << time_step << " s" << endl; 

  /* calculate temperature in each cell */
  F.computeTemp( time );
  double total_loading(0.0);
  
  // Look at the decay of heat over time
  while ( time < 20*day )
    {
      /* Define Creep rate and load fault */
      CreepEvent this_creep( vp , slip_deficit_rate , time_step , time , F );

      /* Add record of the creep to the fault */
      F.add_creep_record( this_creep.getCreepslips_ptr()  , time );

      total_loading += vp*time_step;
      time += time_step;
      it++;

      F.computeStress();
      F.computeTemp( time  );
      write_temp( it, string(argv[0]), time, F ); // write temp
      write_stress( it, string(argv[0]), time, F );

      time_step = dtmax;
      time_step = F.getSlipDefRate( slip_deficit_rate, vp, 
				    time_step, maxcreepslip, 
				    time , true );
      cout << "dt = " << time_step << " s" << endl; 
      write_sliprate( it, string(argv[0]), time, F , slip_deficit_rate , vp ); // write slip rate
      ofile << it << " " << time/hour << endl;
   
      double avgT = F.getAvgTemp();
      cout << "Time = " << time/hour 
	   << " hrs, Avg. T = " << avgT -273 
	   << " ^o C, Max T = " << F.maxTemp()-273 << "o C\n";
      

    }

  write_slipPart( "txtfiles/final_slippart.txt" , F, total_loading );

  ofile.close();
}
//======================================================================
//######################################################################
//######################################################################
