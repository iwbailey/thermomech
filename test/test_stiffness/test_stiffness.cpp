/**
 * Output the change in stress due to slip in a single cell
 */
#include "units.h"
#include "stiffness_matrix.h"
#include <iostream>
#include <stdlib.h> // abs

using namespace std;

//######################################################################
int main(int argc, char * argv[]) {

  cerr << "Program: " << argv[0] << " No. Args: " << argc << endl;

  /* Set the required variables */
  int nx(128), nz(32); // dimensions of the fault grid
  double dx = (70.0*Units::km)/nx, dz = (17.5*Units::km)/nz; // dimensions of a slip surface
  //  double z0 = 0.5*dz;
  double z0 = 5.0*Units::km;
  double rigidity = 30*Units::GPa; // halfspace rigidity
  int iSlip(64), jSlip(8); // location of the slipping cell
  double du(1.0*Units::m); // amount of slip

  cerr << "dx = " << dx << endl;
  cerr << "nx = " << nx << endl;
  cerr << "dz = " << dz << endl;
  cerr << "nz = " << nz << endl;
  cerr << "mu = " << rigidity << endl;
  cerr << "iSlip = " << iSlip << endl;
  cerr << "jSlip = " << jSlip << endl;
  cerr << "du = " << du << endl;

  /* Initialize the stiffness matrix */
  StiffnessMatrix F( dx, nx, z0, dz, nz, rigidity );

  /** Change the output format */
  cout.setf(ios::scientific);
  cout.precision(15);

  /** Loop through all cells on the fault and compute stress change) */
  for( int i = 0; i<nx; i++ )
    for (int j = 0; j < nz; j++)
      {
	/* Along strike difference */
	int ik = abs( i-iSlip );

	/* Get stress change */
	double dtau =  du*F(ik, jSlip, j);

	/* Output the cell center coordinates and stress change */
	cout << ((double)i + 0.5)*dx / Units::km
	     << "\t " << ((double)j + 0.5)*dz / Units::km
	     << "\t " << dtau / Units::MPa
	     << endl;

        if( ik==0 && jSlip - j == 0 )
          cerr << "dtau in slipping cell: " << dtau/Units::MPa << " MPa\n";
  }
}
//======================================================================
//######################################################################
//######################################################################
