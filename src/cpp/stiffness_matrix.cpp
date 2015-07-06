#include "stiffness_matrix.h"


//------------------------------------------------------------------------------
/**
 * Initialize the stiffness matrix using strn().
 * nl, nd are global parameters defined previously
 */
void StiffnessMatrix::stiffness_init(
				     const double &dx, // along-strike width of a slip cell
				     const double &dz, // vertical dimension of a slip cell
                                     const double &z0, // depth at the center of the top cell
				     const double &mu ) // rigidity
{
  double xd2 = 0.5*dz; // cell half-depth
  double xw2 = 0.5*dx; // cell half-width

  /* Loop through all depth cells on the source */
  for( int j1=0; j1<_nz; j1++ )
  {
    /* Depth of source */
    double z1 = z0 + ( (double)j1+0.5 )*dz;

    /*Loop through depth cells of the receiver */
    for( int j2=0; j2<_nz; j2++ )
      {
	/*Depth of receiver*/
	double z2 = z0 + ( (double)j2+0.5 )*dz;

	/* Difference in strike index btw source and receiver. */
	for (int i=0; i<_nx; i++ )
	  {
	    double diffx = (double)i * dx;
	    set_stiff( i, j1, j2, 2 * mu * strn(xw2, xd2, z1, diffx, z2) );
	  }
      }
  }
}
//------------------------------------------------------------------------------
/**
 * Calculates the shear strain for a given point source due to slip in a cell
 * at depth x3.
 */
double StiffnessMatrix::strn(const double &xw2, // Half width of a cell
			     const double &xd2, // Half height of a cell
			     const double &x3,  // Depth of slipping cell center
			     const double &y1,  // Distance along strike for stress calculation
			     const double &y3 )// depth at which to calculate stress
{
  /*
    Copied from Ben-Zion Fortran code
    Ref:
    Chinnery, bssa, 1963 y1 par flt y2 nor flt y3 dwn analytic soln
    pers. com. 1983 xlam=xmu x1=x2=y2=0 tensor strain e12 on fault
    plane x1(=0),x3 (source cell center); y1,y3 (stress calculation
    point)
  */
  double sgn[5], fx1[5], fx3[5];
  double cpi=1.0/(4.0*3.1415926);
  double c23=2.0/3.0;

  sgn[1] = 1.0; // Fortran convention, don't use zero location
  sgn[2] = -1.0;
  sgn[3] = -1.0;
  sgn[4] = 1.0;

  fx1[1]= +xw2; // xw2,xd2 are half width, half height
  fx1[2]= +xw2;
  fx1[3]= -xw2;
  fx1[4]= -xw2;
  fx3[1]= x3 + xd2;
  fx3[2]= x3 - xd2;
  fx3[3]= x3 + xd2;
  fx3[4]= x3 - xd2;
  double e12 = 0.0; // initialize zero shear strain

  for (int i=1; i <= 4; ++i)
    { // label 20
      double t = fx1[i]-y1;
      double q = fx3[i]-y3;
      double p = fx3[i]+y3;
      double tt = t*t, qq = q*q, pp = p*p;
      double s1 = sqrt(tt+qq);
      double s2 = sqrt(tt+pp);
      double s1q = s1 + q;
      double s2p = s2 + p;
      double f1 = 1./(s1*s1q) + 1./(s2*s2p);

      double f2 = (s2/4.0 + q)/(s2*s2p*s2p) - (pp-qq)*(2.0*s2 + p)
	/(2.0*s2*s2*s2*s2p*s2p);
      double f4 = q/(s1*(s1+t)) + p/(s2*(s2+t));
      double etmp = c23*t*(f1 + f2) + f4/2.0;

      e12 += sgn[i]*cpi*etmp;
    }
  // ----- return the calculated shear strain.
  return e12;
}




