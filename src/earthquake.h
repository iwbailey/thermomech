#ifndef _EARTHQUAKE_H_
#define _EARTHQUAKE_H_

/**
 * Functions relating to the ocurrence of an earthquake in the model
 */

#include <vector>
#include <math.h> // for log10 and sqrt
#include "units.h"

#include <iostream> // TMP Debugging

//------------------------------------------------------------------------------
class Earthquake
{

public:

  /**
     Constructor takes all the parameters, does no computation on the fault
  */
  Earthquake(
             const int &hypo_i , const int &hypo_j ,// coordinates of the hypocenter
             const double &time , // time now
             const double &cellArea,
             const std::vector<double> &slipDef0,
             const std::vector<double> &slipDef1,
             const std::vector<double> &stress0,
             const std::vector<double> &stress1 ):
    _t(time),
    _hi(hypo_i),
    _hj(hypo_j),
    _cellArea(cellArea)
  {

    unsigned int nCells = slipDef0.size();
    /* Record only for cells where the slip deficity changed, i.e., the rupture*/
    for( unsigned int k=0; k<nCells; k++ )
      if( slipDef1[k] !=  slipDef0[k] )
        {
          _cellIdxs.push_back( k );
          _cellSlips.push_back( slipDef0[k] - slipDef1[k] );
          _cellStressDrops.push_back( stress0[k] - stress1[k] );
        }
  }

  double failureTime(){ return _t; }// time
  double ruptureArea(){ return _cellArea*_cellIdxs.size(); } // Rupture area
  double potency(); // return the scalar potency
  double potMagnitude(); // Magnitude based on potency
  double moment( const double &rigidity );
  double momMagnitude( const double &rigidity );
  double staticStressDrop(); // static stress drop for the earthquake

private:
  double _t; // The time of earthquake rupture
  int _hi , _hj; // The hypocenter coordinates
  double _cellArea; // Area of a slip cell
  std::vector<unsigned int> _cellIdxs;
  std::vector<double> _cellSlips;
  std::vector<double> _cellStressDrops;

};


//==============================================================================
#endif /* EARTHQUAKE_H_ */
