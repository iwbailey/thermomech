#ifndef ROMBERG_INTEGRATION_H_
#define ROMBERG_INTEGRATION_H_

/**
   Numerical recipes objects/functions used for numerical integration
*/
#include <vector>
#include <iostream>
#include <math.h>

//------------------------------------------------------------------------------
/**
    Polynomial interpolation object
*/
struct Poly_interp {

  const double *xx, *yy; // pointers to arrays defining function evaluation
  int n; // number of points in x and y
  int mm; // Number of points to use in the interpolation
  double dy; // error estimate of last interpolation

  /**
     Constructor
  */
  Poly_interp(
	      std::vector<double> &x, // points of evaluation
	      std::vector<double> &y, // evaluated values of function
	      int m ) // number of points to use in the interpolation
  : n(x.size()), mm(m), xx(&x[0]), yy(&y[0]), dy(0.0) {}

  /**
      Return interpolated value of f(x) at x by mm-point interpolation
      using values of xx and yy from index j1 to j1+mm-1
  */
  double rawinterp( int j1, double x );
};

//------------------------------------------------------------------------------
/**
    Routine implementing extended trapezoidal rule
*/
struct Trapzd {

  int n; // resolution of integration, number of times to divide by 2
  double a, b; // limits of integration
  double s; // current value of integral
  double (*func)(double); // function to integrate, must evaluate a double

  /** Default constructor */
  Trapzd() {}

  /** Constructor takes a function and limits of integration */
  Trapzd( double (*funcc)(double), const double aa, const double bb ):
    func(funcc), a(aa), b(bb), n(0) {}

  /** nth stage of refinement of the extended trapezoidal rule */
  double next();

};

//------------------------------------------------------------------------------
struct Romberg{

  /* Static members, should be defined at the beginning of main cpp
     file. Good initial values are jmax=20, k=5, eps = 1.0e-10 */
  static int JMAX; // limits the total number of steps
  static int K; // number of points used in the interpolation, k=2 is Simpson's rule
  static double EPS; // accuracy of the integral as a proportion of its value

  /**
      Integral of the function from a to b, using Romberg's method of
      order 2K
  */
  double integrate( double (*func)(double), double a, double b );

};

//------------------------------------------------------------------------------

#endif /* ROMBERG_INTEGRATION_H_ */
