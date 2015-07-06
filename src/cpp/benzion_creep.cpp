#include "benzion_creep.h"

//------------------------------------------------------------------------------
/**
   Return the temperature independent creep parameter in Ben-Zion
   (1996), JGR, eqn 2

   x = distance from edge of fault (note difference from BZ1996 where it is along strike coord)
   vPl = plate velocity
   strengthz = background strength at that depth
   xBD = creeping boundary width
*/
double creep_bz96( const double &x, const double &vPl, const double &strengthz,
		   const double &xBD )
{
  double n = 3.0; // stress exponent

  /*Stress ratio of 4 means stress required for vPl at x=0 is 0.25
    strength at xBD*/
  double ln_tauRatio = log(4.0);

  return ( vPl/pow( strengthz, n) ) * exp( n*ln_tauRatio*(1 - x/xBD) );
}

//------------------------------------------------------------------------------
void faultcreep_bz96(
                     std::vector<double> &creepcoeffsBZ96,
                     const int &nL,
                     const int &nD,
                     const double &cellLength,
                     const double &cellHeight,
                     const double &xBD,
                     const double &zBD,
                     const double &strainRate,
                     const double &strength_zBD )
{
  /* Initialize the container for the fault */
  creepcoeffsBZ96.resize( nL*nD );

  /* Calculate the fault depth */
  double faultDepth = (double)(nD)*cellHeight;

  /* Loop along strike */
  for( int i=0; i<nL; i++ )
    {
      /* Compute the distance from the side */
      int nToSide = std::min( i, nL-i-1 );
      double sideDist = ((double)nToSide + 0.5)*cellLength;

      /* Loop down dip */
      for( int j=0; j<nD; j++ )
        {
          /* Distance from the bottom */
          double bottomDist = faultDepth - ((double)j + 0.5)*cellHeight;

          /* Keep the side or bottom result, whichever is lowest strength */
          creepcoeffsBZ96[i*nD+j] = std::max( creep_bz96( sideDist, strainRate,
                                                          strength_zBD, xBD ),
                                              creep_bz96( bottomDist, strainRate,
                                                          strength_zBD,
                                                          faultDepth-zBD ) );
        } // - END depths loop
    } // - END strike loop
  return;
}


// benzion_creep.cpp ends here
