/**
    Source file for the member functions of FaultCoolingModel
*/
#include "fault_cooling.h"

//------------------------------------------------------------------------------
/* Imediate temperature from heat rate over a certain area */
double FaultCoolingModel::temperature( const double &heatRate,
				       const double &area )
{
  return 0.5*heatRate/( area*faultWidth*density*specHeat );
}

//------------------------------------------------------------------------------
/* Temperature from heat rate */
double FaultCoolingModel::temperature( const double &heatRate, const double &area,
				       const double &t1, const double &t2,
				       const double &tNow )
{
  return temperature(heatRate, area)*integTempDecay( t1, t2, tNow );
}

//------------------------------------------------------------------------------
/* Temperature decay function for generated heat  */
double FaultCoolingModel::integTempDecay(
		      const double &tSlipStart, // start time of heat generation
		      const double &tSlipEnd, // end time of heat generation
		      const double &tNow ) // measurement time
{
  /* Check what time to return heat for */
  if( tNow <= tSlipStart){
    /*Before t=0*/
    return 0.0;
  }
  else{
    /* After t=0 */
    double tSinceStart = tNow - tSlipStart;
    double tSinceEnd = fmax(0.0, tNow - tSlipEnd);

    if( fabs( tSinceEnd-tSinceStart ) < Romberg::EPS*tSinceEnd ){
      /* Slip pulse */
      return (tSinceStart-tSinceEnd)*tempDecay( 0.5*(tSinceStart+tSinceEnd) );
    }
    else{
      /* Period of slip */
      return R.integrate( &(tempDecay), tSinceEnd, tSinceStart );
    }
  }
}

//------------------------------------------------------------------------------
/* Temperature decay function */
double FaultCoolingModel::tempDecay( double tSinceSlip )
{
  /* Check greater than 0.0 */
  if( fabs(tSinceSlip) >= Romberg::EPS ){
    double tmp = 0.25*faultWidth/sqrt( diffusivity*tSinceSlip );
    return erf(tmp) - erf(-tmp);
  }
  else{
    /* Case where tSinceSlip = 0.0: erf(Inf) - erf(Inf) = 1 - -1 = 2 */
    return 2.0;
  }
}
//------------------------------------------------------------------------------
