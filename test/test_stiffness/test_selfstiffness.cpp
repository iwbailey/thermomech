/**
 * Output the self stiffness column to stdout
 *
 *
 */
#include "units.h"
#include "stiffness_matrix.h"
#include <iostream>

using namespace std;

//######################################################################
int main(int argc, char * argv[]) {

  cerr << "Program: " << argv[0] << " No. Args: " << argc << endl;

  /* Set the required variables */
  int nx(128), nz(32); // dimensions of the fault grid
  double dx = (70.0*Units::km)/nx, dz = (17.5*Units::km)/nz; // dimensions of a slip surface
  double rigidity = 30*Units::GPa; // halfspace rigidity

  cerr << "dx = " << dx << endl;
  cerr << "nx = " << nx << endl;
  cerr << "dz = " << dz << endl;
  cerr << "nz = " << nz << endl;
  cerr << "mu = " << rigidity << endl;

  /* Initialize the stiffness matrix */
  StiffnessMatrix F( dx, nx, dz, nz, rigidity );

  /** Change the output format */
  cout.setf(ios::scientific);
  cout.precision(15);

  /* Output only the self stiffness  */
  for (int j = 0; j<nz; j++)
    cout << ((double)j+0.5)*dz/Units::km << "\t " << F(0, j, j) << endl;

}
//======================================================================
//######################################################################
//######################################################################
