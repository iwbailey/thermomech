#include "side_creep.h"

//------------------------------------------------------------------------------
/**
   Adjust Arrhenius amplitude and activation energy within a distance xBD from
   the side of the fault so that...
   1. there is no temperature dependent creep for x<xBD and x>(L-xBD)
   2. the strength profile as a function of x reflects the strength profile as a
      function of (W-z) at x=xBD, where W is the fault width, and z is the depth.
*/
void sidecreepmask( std::vector<double> &A,
                    std::vector<double> &E,
                    const std::vector<double> &n,
                    const std::vector<double> &T,
		    const int &nL, const int &nD,
                    const double &cellLength,
                    const double &xBD )
{
  /* Check dimensions of input arguments */
  if( A.size() != nL*nD || E.size() != nL*nD ){
    std::cerr << "ERROR: size of A and E must be nL*nD" << std::endl;
    std::terminate();
  }

  int iSide = 0;
  while( ((double)iSide + 0.5)*cellLength < xBD ) iSide++;

  /* Loop through along strike grid coordinates each grid point */
  for( int i=0; i<iSide; i++ )
    {
      //double sideDist = ((double)i + 0.5)*cellLength;
      int iRev = nL-i-1;

      /* Get the reference cell to mimic */
      int kRef = iSide*nD + nD-i-1;
      CreepLaw C1left( A[kRef], E[kRef], n[kRef] );
      double A1left = C1left.temperatureTerm( T[kRef] );

      /* Repeat for the right hand side */
      kRef = (nL-1-iSide)*nD + nD-i-1;
      CreepLaw C1right( A[kRef], E[kRef], n[kRef] );
      double A1right = C1right.temperatureTerm( T[kRef] );

      /* Loop through down-dip grid coordinates */
      for( int j=0; j<nD; j++ )
        {
          /* Set the alternative A to mimic what is there */
          CreepLaw C2left( A[i*nD+j], E[i*nD+j], n[kRef] );
          double A2left = C2left.temperatureTerm( T[i*nD+j] );

          CreepLaw C2right( A[iRev*nD+j], E[iRev*nD+j], n[kRef] );
          double A2right = C2right.temperatureTerm( T[iRev*nD+j] );

          /* Turn off temperature dependence by setting E=0 */
          E[i*nD+j] = 0.0;
          E[iRev*nD+j] = 0.0;

          /* Choose whichever is weaker */
          A[i*nD+j] = std::max( A1left, A2left );
          A[iRev*nD+j] = std::max( A1right, A2right );

        } // END Loop through depth

    } // END Loop through bottom

  return;
}

// side_creep.cpp ends here
