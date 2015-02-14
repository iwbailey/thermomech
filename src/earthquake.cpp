#include "earthquake.h"

//------------------------------------------------------------------------------
/** Compute the scalar potency */
double Earthquake::potency()
{

  double potency = 0.0;
  for( unsigned int k=0; k<_cellIdxs.size(); k++ )
    potency += _cellSlips[k]*_cellArea;
  return potency;
}

//------------------------------------------------------------------------------
/** Compute the magnitude based on the potency scaling relation */
double Earthquake::potMagnitude()
{
  /* Use the empirical relation of Ben-Zion and Zhu (GJI, 2002), quadratic form */

  /* Quadratic terms */
  double a = 0.0612;
  double b = 0.988;

  /* Make sure we have the correct units for the scaling relation*/
  double c = -4.87 - log10( potency() /(Units::km*Units::km*Units::cm));

  /* Use the quadratic formula */
  double mag = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);
  return mag;
}

//------------------------------------------------------------------------------
/** Compute the scalar moment */
double Earthquake::moment( const double &rigidity )
{
  return rigidity*potency();
}

//------------------------------------------------------------------------------
/** Compute the magnitude based on the Hanks & Kanamori scaling relation */
double Earthquake::momMagnitude( const double &rigidity )
{
  /* log10 M0 = 1.5 mW + 16.1 */
  double logM = log10( moment( rigidity ) / (Units::N*Units::m) );
  return (logM - 9.1)/1.5;
}

//------------------------------------------------------------------------------
/** Compute the static stress drop, i.e. the average drop in stress for all
    ruptured cells*/
double Earthquake::staticStressDrop()
{
    unsigned int nFailed = _cellIdxs.size();

    /* Add up total stress drops */
    double stressdrop =0.0;
    for( unsigned int i=0; i< nFailed; ++i )
      stressdrop += _cellStressDrops[i];

    /* Return the average */
    return stressdrop / nFailed;
  }

