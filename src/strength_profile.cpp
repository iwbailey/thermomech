#include "strength_profile.h"

//------------------------------------------------------------------------------
/* Define profile on a fault grid, column ordered */
std::vector<double>
StaticStrengthProfile::faultgrid_co(
                                    const int &nL, const int &nD,
                                    const double &cellHeight
                                    ) const
{
  /* Set up the output containter */
  std::vector<double> faultStrengths( nL*nD );

  /*Loop through all depths */
  for( int j=0; j<nD; j++ )
    {
      /* Calculate the depth */
      double depth = ((double)j + 0.5)*cellHeight;

      /* Fill in the temerature for all cells in this depth row */
      for( int i=0; i<nL; i++ ) faultStrengths[i*nD+j] = (*this)(depth);
    }

  /* Return the vector */
  return faultStrengths;
}
//------------------------------------------------------------------------------

