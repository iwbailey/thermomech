/**
   Code for functions associated with the class ThermoFault in thermofault.h
*/
#include "thermofault.h"

//------------------------------------------------------------------------------
/**
    Get a vector array of temperature at all points on the fault
*/
void ThermoFault::getTemperature(
			   std::vector<double> &temperature, // OUT column ordered array of temperature on the fault
			   const double &timeNow
			   )
  const
{
  /* Start with the background temperature */
  temperature = _initTemperature;

  /* Loop through earthquakes stored for the fault */
  for( unsigned int iEqk = 0; iEqk < _heatHistory.size(); iEqk ++ ){

    /* Compute the heat decay part of the cooling equation */
    double heatDecay = _FCM->integTempDecay( (_heatHistory[iEqk]).t1,
                                            (_heatHistory[iEqk]).t2, timeNow );

    /* Add on the individual heat rate for all cells involved in the
       earthquake */
    for( unsigned int iCell=0; iCell < _heatHistory[iEqk].failedCellIdxs.size(); iCell++ ){
      int k = _heatHistory[iEqk].failedCellIdxs[iCell];
      temperature[k] += _FCM->temperature( _heatHistory[iEqk].cellHeatRates[iCell],
                                          _cellHeight*_cellLength ) * heatDecay;
    }
  }

  return;
}


//------------------------------------------------------------------------------
/**
    Based on the total slip for each cell record the heat generated in an
    earthquake */
void ThermoFault::recordHeatEvent( const std::vector<double> &heatFromSlip,
                             const double &t1, const double &t2 )
{
  /* Containers for which cells failed and what their heat rates were */
  std::vector<int> failedCellIdxs;
  std::vector<double> heatRates;

  /* Record the heat rate for each failed cell */
  for( int k=0; k<_nL*_nD; k++ )
    {
      /* Check if there was any slip in this cell */
      if( heatFromSlip[k] > 0.0 ){
        failedCellIdxs.push_back( k );

        /* Distributed heat equally over the time of slip */
        heatRates.push_back( heatFromSlip[k]/ (t2-t1) );
      }
    }

  /* Add heat event to the containter for the fault used in temperature calc*/
  _heatHistory.push_back( HeatEvent( t1, t2, failedCellIdxs, heatRates ) );
}

//------------------------------------------------------------------------------
/**
 * Compute the cascading failure of cells given an initial "hypocenter" cell
 * already in a critical state.  Update the slip deficit, stress and time.
 * Output the distribution of total slip on the fault.
 *
 IN/OUT
  stress : current stress on fault, updated by function
  time : current time, updated according to the max amount of slip and the slip
  velocity
 IN
  iHypo : along strike index of the hypocenter cell
  jHypo : down-dip index of the hypocenter cell
  slipVelocity: constant slip velocity of seismic slip in the model
 OUT
  totalSlip : slip in each cell of the fault during the earthquake

*/
void ThermoFault::computeEarthquake( std::vector<double> & stress,
                                     double &time,
                                     const int &iHypo, const int &jHypo,
                                     const double &slipVelocity,
                                     std::vector<double> &totalSlip,
                                     std::vector<double> &totalHeat )
{
  /* Check that the hypocenter is critical */
  int kHypo = cellidx(iHypo,jHypo);
  if( stress[kHypo] < _strengthStatic[kHypo] )
    throw std::runtime_error("Hypocenter is not critical");

  /* Container for recording the total slip */
  totalSlip.resize( _nL*_nD );
  std::fill( totalSlip.begin(), totalSlip.end(), 0.0 );

  /* Container for recording the total heat*/
  totalHeat.resize( _nL*_nD );
  std::fill( totalHeat.begin(), totalHeat.end(), 0.0 );

  /* Set up the first critical cell */
  std::vector<bool> isStatic( _nL*_nD, true );
  std::vector<int> criticalIdxs( 1, kHypo );
  int nCritical = 1;

  /* Compute stress redistribution and slip until all cells are below
     strength */
  while( nCritical > 0 ){

    /* Loop through failed cells */
    for( unsigned int iCrit=0; iCrit < nCritical; iCrit++ )
      {
        /* Get the cell index*/
        int k = criticalIdxs[iCrit];

        /* Compute the stress drop and the amount of slip*/
        double stressdrop = stress[k] - cellArrestStress(k);
        double slip = -1*stressdrop / _Stiffness( 0, k%_nD, k%_nD);

        /* Compute the amount of heat, assuming linear decrease in stress */
        double slidingStress = cellArrestStress(k) + 0.5*stressdrop;
        double heat = _cellHeight*_cellLength * slip * slidingStress;

        /* Update the containers for this cell */
        _slipDeficit[k] -= slip;
        totalSlip[k] += slip;
        totalHeat[k] += heat;
        isStatic[k] = false;
      }

    /* Recompute the stress */
    getStress( stress );

    /* Find any failed cells */
    nCritical = findCriticalCells( criticalIdxs, stress, isStatic );
  }

  /* Calculate the time taken based on the maximum slip amount*/
  double tEnd = time +
    *(std::max_element( totalSlip.begin(), totalSlip.end() ))/slipVelocity;

  /* Store a record of heat generated for each cell that slipped */
  recordHeatEvent( totalHeat, time, tEnd );

  /* Update time */
  time = tEnd;
}
