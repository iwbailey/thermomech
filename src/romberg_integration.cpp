/* Source code for definitions in romberg_integration.h */
#include "romberg_integration.h"

//--------------------------------------------------------------------------------
/**
    Interpolated value y at x by mm point interpolation
*/
double Poly_interp::rawinterp( int j1, double x )
{
  int i;
  double y, den, ho, hp, w;
  std::vector <double> c(mm), d(mm);

  const double *xa = &xx[j1], *ya = &yy[j1];

  /* Find index of nearest table entry between j1 and j1+mm*/
  int iS = 0;
  double dif;
  double minDif = fabs( x-xa[0] );
  for( i=0; i<mm; i++)
    {
      dif = fabs(x-xa[i]);
      if( dif < minDif ){
	iS = i;
	minDif = dif;
      }
      c[i] = ya[i];
      d[i] = ya[i];
    }

  /* Initial approximation to y */
  y = ya[iS--];

  /* Loop over current c's and d's and update them */
  int m;
  double dxLow, dxHigh;
  for( m=1; m<mm; m++ )
    {
      for( i=0; i<mm-m; i++ )
	{
	  dxLow = xa[i] - x;
	  dxHigh = xa[i+m] - x;
	  w = c[i+1] - d[i];

	  /* Error if two xa's are identical */
	  den = dxLow - dxHigh;
	  if( den == 0.0 ){
	    std::cerr << "poly_interp error\n";
	    std::terminate();
	  }
	  den = w/den;
	  d[i] = dxHigh*den;
	  c[i] = dxLow*den;
	}
      if( 2*iS < mm-m ) dy = c[iS+1]; // forking up
      else dy = d[iS--]; // forking down, update iS
      y += dy;
    }
  return y;
}

//------------------------------------------------------------------------------
/**
   Return the value of the integration at an increased resolution
   level to the currently stored value
 */
double Trapzd::next()
{
    double x, tnm, sum, del;
    int it, j;

    /* Increase the resolution of integration */
    n++;
    if( n==1 ){
      /* Simple averaging */
      s = (b-a)*0.5*(func(a)+func(b));
      return s;
    }
    else{
      /* Get number of intervals */
      for( it=1, j=1; j<n-1; j++ ) it <<= 1; // it = it*2^1
      tnm = (double) it;

      /* New spacing of points to be added */
      del = (b-a)/tnm;

      /* Advance x to first point */
      x = a + 0.5*del;

      /* Sum f(x) at all x points */
      for( sum=0.0, j=0; j<it; j++, x+=del) sum += func(x);

      /* Combine with existing s to get new s */
      s = 0.5*(s+(b-a)*sum/tnm);
      return s;
    }
}
//-----------------------------------------------------------------------------
/** Returns the integral of the function from a to b, using Romberg's
    method of order 2K
*/
double Romberg::integrate( double (*func)(double), double a, double b)
{

  /* Stores for successive tapezoidal approx*/
  std::vector<double> s(JMAX), h(JMAX+1);

  /* Interpolation structure for final integration */
  Poly_interp polint(h,s,K);

  /*x axis of integration, start at 1, decrease by factor 0.25 every
    successive approximation*/
  h[0] = 1.0;

  /* Set up the trapezoidal integration */
  Trapzd t( func, a, b );

  /* Loop through all steps */
  for( int j=1; j<=JMAX; j++ )
    {
      /* Compute next estimate of integral from trapezoidal rule */
      s[j-1] = t.next();

      /* If we have enough points for the interpolation, apply it */
      if( j>= K ){
	/* Interpolate to lim h->0 */
	double ss = polint.rawinterp( j-K, 0.0 );

	/* Check the error */
	if( fabs( polint.dy) <= EPS*fabs(ss)) return ss;
      }

      /* Decrease h */
      h[j] = 0.25*h[j-1];
    }

}

//------------------------------------------------------------------------------
