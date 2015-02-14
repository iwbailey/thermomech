/**
   Code for functions associated with the class Fault in fault.h
*/

#include "fault.h"

//------------------------------------------------------------------------------
/**
   Add plate loading to the slip deficit when creep velocity known
*/
void Fault::loadFault(
                      const double &plateVelocity,
                      const std::vector<double> &creepVel,
                      const double &dtime )
{

  /* Update the slip deficit in all cells */
  for( int k=0; k<(_nL*_nD); k++ )
    {
      /* Compute the amount of slip deficity to add*/
      double slipDeficit = dtime*( plateVelocity - creepVel[k] );

      /* Add the remaining loading */
      _slipDeficit[k] += slipDeficit;
    }

  return;
}

//------------------------------------------------------------------------------
/**
   Calculate the stress based on slip deficit at cell i,j
*/
double Fault::getCellStress(
			    const int &i, // Along strike index
			    const int &j ) // Down dip index
  const
{

  /* Start with initial stress and add contribution from all cells */
  double stress = _initStress[ cellidx(i,j) ];

  /* Loop through all cells on the fault */
  for (int l=0; l< _nD; ++l)
    for (int k=0; k < _nL; ++k)
      {
	/* Get the along strike difference from queried cell */
	int ik = abs(i-k);

	/* Add the influence onto the total stress */
	//  double checked order of indices!!
	stress += _Stiffness(ik,j,l) * ( -1.0*_slipDeficit[ cellidx(k,l) ] );
    }

  /* Check for errors */
  if( isnan(stress) ) throw std::runtime_error("NaN stress");

  return stress;

}

//------------------------------------------------------------------------------
/**
    Get a vector array of stress at all points on the fault
    IN/OUT column ordered array of stress on the fault
*/
void Fault::getStress( std::vector<double> &stress  )
  const
{
  /* Check the dimensions of the vector */
  if( stress.size() != _nL*_nD ){
    stress.resize( _nL*_nD );
  }

  /* Loop through all cells on the fault */
  for (int i=0; i < _nL; ++i)
    for (int j=0; j < _nD; ++j)
      {
	stress[cellidx(i,j)] = getCellStress(i,j);
      }
      // if ( _tau[i*_nd+j] < 0.0 )
      // 	cerr << "WARNING: -ve stress in cell [" << i << "," << j << "]\n";

  return;
}


//------------------------------------------------------------------------------
/**
   Get the aseismic strain rate of slip cell k given stress and temperature
*/
double Fault::getCellStrainRate( const int &k, const double &stress,
				 const double &temperature )
  const
{
  /* Set up the creep law parameters for this cell */
  CreepLaw C( _arrhA[k], _activEnergy[k], _stressExpon[k] );

  /* Get the strain rate at this temp and stress */
  return C.strainRate( stress, temperature );
}

//------------------------------------------------------------------------------
/**
   Get the aseismic slip velocity for slip cell k given stress and temperature
*/
double Fault::getCellCreepVelocity( const int &k, const double &stress,
				    const double &temperature )
  const
{
  return _faultWidth * getCellStrainRate( k, stress, temperature );
}

//------------------------------------------------------------------------------
/**
   Get the aseismic slip velocity for all cells
*/
void Fault::getCreepVelocity( std::vector<double> &creepVel,
                              const std::vector<double> &stress,
                              const std::vector<double> &temperature )
  const
{
  /* Check the in/out vector size */
  if( creepVel.size() != _nL*_nD )
    creepVel.resize(_nL*_nD);

  /* Loop through all cells and get the result */
  for( int k=0; k<_nL*_nD;  k++ )
    creepVel[k] = getCellCreepVelocity( k, stress[k], temperature[k] );

  return;
}

//------------------------------------------------------------------------------
/**
   Estimate time to failure for the entire fault given loading at current rate
   and stress. Note that this does not account for changes to the creep velocity
   caused by changes to the stress. Hence, this function should not be used to
   set the time step unless constrained by an upper bound. */
double Fault::estimateTimeToFailure(
				    const std::vector<double> &stress,
				    const std::vector<double> &creepVel,
				    const double &plateVelocity )
  const
{
  /* Initial time is a long time */
  double minTime = 99e99;

  /* Loop through all cells on the fault */
  for( int i=0; i<_nL; i++ )
    for( int j=0; j<_nD; j++ )
      {
        int iCell = cellidx(i,j);

        /* Amount of stress increase needed to reach failure */
        double stressToFailure = _strengthStatic[iCell] - stress[iCell];

        /* Check if already at failure */
        if( stressToFailure <= 0 ) return 0.0;

        /* Convert the loading rate from all cells on the fault into a stress
           build up rate for this fault, assuming constant loading */
        double stressingRate = 0.0;
        for( int i2=0; i2<_nL; i2++ )
          for( int j2=0; j2<_nD; j2++ )
            {
              double ff = _Stiffness( abs(i-i2), j, j2);
              double loadingRate = plateVelocity - creepVel[ cellidx(i2,j2) ];
              stressingRate += ff*( -1.0*loadingRate );
            }

        /* Compute the time taken to reach failure*/
        double failTime = stressToFailure / stressingRate;

	/* Store the minimum positive rate */
	if( stressingRate > 0 && failTime < minTime ) minTime = failTime;

      }

  /* Result is the minimum failure time */
  return minTime;

}

//------------------------------------------------------------------------------
/**
 * Return number of cells above failure threshold critical
 **/
int Fault::nCriticalCells(
                          unsigned int &iHypo,
                          unsigned int &jHypo,
                          const std::vector<double> &stress )
  const
{
  double maxtau = -99.99;
  int nCrit = 0;

  for( int i=0; i<_nL; i++ )
    for( int j=0; j<_nD; j++ )
      {
        int k = cellidx(i,j);
        double overstress = stress[k] - _strengthStatic[k];

        /* check if critical */
        if( overstress >= 0.0 ){
          nCrit++;

          /* Check if most critical */
          if( overstress > maxtau ){
            maxtau = overstress;
            iHypo = i;
            jHypo = j;
          }
        }
      } // END loop through fault cells

  return nCrit;
}

//------------------------------------------------------------------------------
/**
 * Return which cells are above failure threshold critical
 **/
int Fault::findCriticalCells(
                             std::vector<int> &criticalIdxs,
                             const std::vector<double> &stress,
                             const std::vector<bool> &isStatic )
  const
{
  /* Clear any existing records */
  criticalIdxs.clear();

  /* Record number of critical cells */
  int nCrit(0), k;
  double strength;

  for( int i=0; i<_nL; i++ )
    for( int j=0; j<_nD; j++ ){
      /* Get the index for the vectors */
      int k = cellidx(i,j);

      /* Find static or dynamic strength */
      if( isStatic[k] ){
        strength = _strengthStatic[k];
      }
      else{
        strength = _strengthDynamic[k];
      }

      /* Add the indices of critical cells to the container */
      if( stress[k] >= strength ){
        criticalIdxs.push_back( k );
        nCrit++;
      }
    } // END loop through fault cells

  return nCrit;
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
void Fault::computeEarthquake( std::vector<double> & stress,
                               double &time,
                               const int &iHypo, const int &jHypo,
                               const double &slipVelocity,
                               std::vector<double> &totalSlip )
{
  /* Check that the hypocenter is critical */
  int kHypo = cellidx(iHypo,jHypo);
  if( stress[kHypo] < _strengthStatic[kHypo] )
    throw std::runtime_error("Hypocenter is not critical");

  /* Container for recording the total slip */
  totalSlip.resize( _nL*_nD );
  std::fill( totalSlip.begin(), totalSlip.end(), 0.0 );

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

        /* Update the containers for this cell */
        _slipDeficit[k] -= slip;
        totalSlip[k] += slip;
        isStatic[k] = false;
      }

    /* Recompute the stress */
    getStress( stress );

    /* Find any failed cells */
    nCritical = findCriticalCells( criticalIdxs, stress, isStatic );
  }

  /* Update time based on the maximum slip amount*/
  time += *(std::max_element( totalSlip.begin(), totalSlip.end() ))/slipVelocity;

}

// //------------------------------------------------------------------------------
// /** Compute the heat generated from a period of creep and add to the fault
//     history*/
// void Fault::calcHeatFromCreep( const double &t1, const std::vector<double> &stress1,
//                         const double &t2, const std::vector<double> &stress2,
//                         const std::vector<double> &creepVel,
//                         const double &minSlipRate )
// {
//   /* Container for recording the heat */
//   std::vector<double> heatFromSlip( _nL*_nD, 0.0);

//   /* Loop through all cells */
//   for( int k=0; k<_nL*_nD; k++ ){
//     /* Check if creep rate exceeds lower limit and that creep depends on temp */
//     if( std::fabs( creepVel[k] ) > minSlipRate && _activEnergy[k] != 0.0 ){
//       /* Assume slip at constant rate */
//       double slip = creepVel[k]*(t2-t1);

//       /* Assume stress at average of start and end stress */
//       double stress = 0.5*(stress1[k] + stress2[k]);
//       heatFromSlip[k] = slip*stress*_cellHeight*_cellLength;
//     }
//   }
//   /* Add to the history */
//   recordHeatEvent( heatFromSlip, t1, t2 );
//   return;
// }


//------------------------------------------------------------------------------
// //   /**
// //      Get Max temperature
// //   */
// //   double maxTemp(){
// //     vector<double> :: iterator pos = max_element (
// // 						  _temperature.begin(),
// // 						  _temperature.end()
// // 						   );
// //     return *pos;
// //   }
// //  /**
// //     Get Max background Temperature, used for setting up background stress
// //  */
// //   double maxTbkd(){
// //     double T, maxT(0.0);
// //     for( int i = 0 ; i < _nl*_nd ; i++ ){
// //       T = bkdT(int(i/_nd),i%_nd);
// //       if( T > maxT ) maxT = T;
// //     }
// //     return maxT;
// //   }

// //   /**
// //       Get average temperature
// //   */
// //   double getAvgTemp(){
// //     double sumT = 0;
// //     for( int i = 0 ; i < _nl*_nd ; i++ ) sumT += _temperature[i];
// //     return sumT / (_nl*_nd);
// //   }

// //   /**
// //      Get average stress
// //   */
// //   double getAvgStress(){
// //     double sumS = 0;
// //     for( int i = 0 ; i < _nl*_nd ; i++ ) sumS += _tau[i];
// //     return sumS / (_nl*_nd);
// //   }

// //   /**
// //      Get average slip deficit
// //   */
// //   double getAvgSlipDef(){
// //     double sumSD = 0;
// //     for( int i = 0 ; i < _nl*_nd ; i++ ) sumSD += _slipdef[i];
// //     return sumSD / (_nl*_nd);
// //   }

// //   /**
// //      Get the seismic slip between the time limits t0 and t1
// //   */
// //   double getSeisSlip( const int &il , const int & id,
// // 		      const double &t0, const double &t1 )
// //   {

// //     return _slipHistory[il*_nd+id].getSlip( t0, t1 );
// //   }

// //   /**
// //       Get the total seismic slip
// //   */
// //   double getTotalSeisSlip( const int &il , const int & id ) // coordinates
// //   {
// //     return _slipHistory[il*_nd+id].getSlip( 0.0 ); // get slip since time = 0.0
// //   }

// //   /**
// //      Get the total number of  seismic slip events
// //   */
// //   double getSeisSlipCount( const int &il , const int & id) // coordinates
// //   {
// //     return _slipHistory[il*_nd+id].getSlipCount();
// //   }

// //   /**
// //      Get the creep slip between the time limits t0 and t1
// //   */
// //   double getCreepSlip( const int &il , const int & id,
// // 		      const double &t0, const double &t1 )
// //   {
// //     return _creepHistory[il*_nd+id].getSlip( t0, t1 );
// //   }

// //   /**
// //      Get the total aseismic slip
// //   */
// //   double getTotalCreepSlip( const int &il , const int & id ) // coordinates
// //   {
// //     return _creepHistory[il*_nd+id].getSlip( 0.0 ); // get slip since time = 0.0
// //   }


// //   /**
// //       A cell has failed so failure is now at dynamic levels
// //   */
// //   void become_dynamic( const int &i , const int & j ){ _tauf[i*_nd+j] = taud(i,j); }

// //   /**
// //    * Lock the fault
// //    */
// //   void lock(){ for( int i = 0 ; i < _nl*_nd ; i++ ) _tauf[i] = _taus[i]; }

// //   /**
// //    * Override the temperature - used in routines that check the change
// //    *  in strength due to temp
// //    */
// //   void setTemperature( const int &i , const int & j , const double & T ){
// //     _temperature[i*_nd+j] = T; }

// //   //============================================================
// //   //=======================
// //   /**
// //    * Return true if the locked fault is critical, return indices of
// //    * critical cells
// //    **/
// //   bool is_critical(
// // 		   vector<int> &crit_i ,
// // 		   vector<int> &crit_j ) //  coordinates of critical cells
// //   {
// //     crit_i.clear(); // empty containers
// //     crit_j.clear();
// //     bool isCritical = false;

// //     for( int il = 0 ; il < _nl ; il++ )
// //       for( int id = 0 ; id < _nd ; id++ ){ // loop over fault
// // 	double overstress = _tau[il*_nd+id] - _tauf[il*_nd+id];
// // 	if( overstress >= 0.0 ){ // check if critical
// // 	  isCritical = true;
// // 	  crit_i.push_back( il );
// // 	  crit_j.push_back( id );
// // 	}
// //       }
// //     return isCritical;
// //   }

// //   //============================================================
// //   /**
// //    * Load the entire fault uniformly (except creepmask)
// //    */
// //   void load(
// //             const double &du )// amount of slip to add
// //   {
// //     for( int i = 0 ; i < _nl*_nd ; i++ )
// //       _slipdef[i] += _creepmask[i]*du;
// //   }

// //   /**
// //    * Load an individual cell
// //    **/
// //   void load( const int &i , const int & j , // coordinates
// // 	     double &delu )// INOUT: amount of slip to add
// //   {
// //     delu *= _creepmask[i*_nd+j]; // adjust for creep mask
// //     _slipdef[i*_nd+j] += delu;
// //   }

// //   /**
// //      load a fault given loading rate and time step
// //   */
// //   void load( vector<double> &loadrate , const double & dt )
// //   {
// //     for( int i = 0 ; i < _nl*_nd ; i++ )
// //       _slipdef[i] += _creepmask[i]*loadrate[i]*dt;
// //   }


// //   //============================================================
// //   /**
// //    *  Unload a cell
// //    **/
// //   void add_slip( const int &i , const int & j , // coordinates
// // 		 const double &delu ) // amount of slip deficit to remove
// //   { _slipdef[i*_nd+j] -= delu; }

// //   //============================================================
// //   /**
// //    * Overloaded copy of above Compute the stress from slip deficit and
// //    *	return an array with stress values
// //    **/
// //   void computeStress( vector<double> &stress )
// //   {
// //     stress.clear();
// //     for (int i=0; i < _nl; ++i) // obs cell length index
// //       for (int j=0; j < _nd; ++j){ // obs cell depth index

// // 	//_tau[i*_nd+j] = 0.0; // compute entirely from slip deficit
// // 	_tau[i*_nd+j] = _tau_bkd[i*_nd+j]; // start with background level

// // 	for (int l=0; l< _nd; ++l) // source cell depth index
// // 	  for (int k=0; k < _nl; ++k){ // source cell length index

// // 	    int ik = abs(int(i)- int(k)); //--- along strike difference
// // 	    _tau[i*_nd+j] += ff(ik,j,l) * ( -_slipdef[k*_nd+l] ); //  double checked order of indices!!
// // 	  } // end k and l
// // 	stress.push_back(_tau[i*_nd+j]);

// // 	if ( _tau[i*_nd+j] < 0.0 )
// // 	  cerr << "WARNING: -ve stress in cell [" << i << "," << j << "]\n";

// //       }// end i and j
// //     return;
// //   }


// //   /**
// //    * Compute the stress from slip deficit
// //    **/
// //   void computeStress()
// //   {
// //     vector<double> tmpstress;
// //     computeStress( tmpstress );
// //     return;
// //   }

// //   //============================================================

// //   /**
// //    * Get the temperature at (i,j) which is sum of background temp and
// //    * temp from heat generation on the fault Return average fault
// //    * temperature
// //    **/
// //   double computeTemp(
// // 		     const double & time  // the time now
// // 		     )
// //   {

// //     double oldT = getAvgTemp();

// //     for (int il=0; il < _nl; ++il)
// //       for (int id=0; id < _nd; ++id)
// // 	{
// // 	  /* 2. add on the heat from the last q eqk events evaluated
// // 	     at the current time */
// // 	  double eqkT = _creepmask[il*_nd+id] *
// // 	    _slipHistory[il*_nd+id].getTempFromSlip(
// // 						    time - _eqktcut ,
// // 						    time ,
// // 						    _cellLength ,
// // 						    _cellDepth ,
// // 						    _thick ,
// // 						    taud(il,id) ,
// // 						    _h );


// // 	  /* 3. add on the heat from the past creep */
// // 	  double creepT = _creepmask[il*_nd+id]*
// // 	    _creepHistory[il*_nd+id].getTempFromSlip(
// // 						     time - _creeptcut ,
// // 						     time ,
// // 						     _cellLength ,
// // 						     _cellDepth ,
// // 						     _thick ,
// // 						     _h );

// // 	  double Tij = bkdT(il,id) + creepT + eqkT;

// // 	  //if( Tij > bkdT(il,id) )
// // 	  //  cout << il << " " << id << " "
// // 	  //	 << bkdT(il,id) << " " << creepT << " " << eqkT << endl;

// // 	  if( Tij < bkdT( il,id ) )
// // 	      cerr << "WARNING: Negative temp from slip at " << il << "," << id << endl;



// // 	  _temperature[il*_nd+id] = Tij;
// // 	}

// //     //cerr << "Avg T change = " <<
// //     return getAvgTemp() - oldT;
// //   }

// //   //============================================================

// //   double computeActEnergy(
// // 			  vector<double> &AE0  // Normal activation energy
// // 			  , vector<double> &AE1  // Activation energy after eqk
// // 			  , const double &t_h // healing time
// // 			  , const double &t_now // time now
// // 			   )
// //   {
// //     /**
// //        Compute activation energy based on time since last earthquake
// //      **/

// //     /* Loop through all cells */
// //     for( int il=0 ; il < _nl ; il++ )
// //       for( int id=0 ; id < _nd ; id++ )
// // 	{
// // 	  /* Compute time since last eqk */
// // 	  double teq = _slipHistory[il*_nd+id].getLastSlipTime();

// // 	  if( teq < 0 || t_now - teq > t_h )
// // 	    /* Case where earthquake has no effect */
// // 	    _actEnergy[il*_nd+id ] = AE0[il*_nd+id];
// // 	  else if( teq == t_now  )
// // 	    _actEnergy[il*_nd+id ] = AE1[il*_nd+id];
// // 	  else{
// // 	    _actEnergy[il*_nd+id] = AE1[il*_nd+id] +
// // 	      (t_now - teq) * (AE0[il*_nd+id] - AE1[il*_nd+id])/t_h;
// // 	  }
// // 	}
// //     return 0.0;
// //   }

// //   //============================================================
// //   /**
// //    * Distribute stress from change slip deficit
// //    **/
// //   void distributeStress( const int & ilc , const int &idc , const double delslip )
// //   {
// //     for (int il=0; il<_nl; ++il)
// //       for (int id=0; id<_nd; ++id)
// // 	{
// // 	  int ik = abs(il - ilc);
// // 	  _tau[il*_nd+id] += ff(ik,id,idc) * delslip;
// // 	}
// //   }

// //   //============================================================
// //   /**
// //    * Get creep vel based on current stress and temp
// //    **/
// //   void getInstSlipRate(
// // 		       vector<double> &vc // OUT: column ordered vector
// // 		       )
// //   {
// //     /* make sure container is right size */
// //     vc.resize( _nl*_nd );

// //     /* Compute the slip rate and time step */
// //     for (int i = 0; i <_nl*_nd; ++i)
// //       {
// //         /* compute creep velocity */
// // 	vc[i] = creepvelocity( // see diffusion.h
// // 			       _thick,
// // 			       _tau[i] , // the stress field
// // 			       _temperature[i] ,// the temperature field
// // 			       _actEnergy[i] ,
// // 			       _h->arrhA() ,
// // 			       R_g );

// //       }

// //     return ;
// //   }


// //   //============================================================
// //   /**
// //    * Get instantaneous creep deformation rate based on
// //    *	temperature and current stress, return time step
// //    **/
// //   double getInstSlipDefRate(
// // 			    vector<double> &slipdef_rate // OUT: column ordered vector
// // 			    , const double &vp  // plate velocity
// // 			    , const double &dtmax // max time step
// // 			    , const double &umax //max allowed creep
// // 			    , bool isTemp // flag = true to use current temp, otherwise background
// // 			    )
// //   {
// //     double dt; // timestep
// //     double maxVc( 0.0 );

// //     /* make sure container is right size */
// //     slipdef_rate.resize( _nl*_nd );

// //     /* Compute the slip rate and time step */
// //     for (int i = 0; i <_nl*_nd; ++i)
// //       {
// //         /* compute creep velocity */
// // 	double vc;

// // 	if( isTemp ){
// // 	  vc = creepvelocity( // see diffusion.h
// // 			     _thick,
// // 			     _tau[i] , // the stress field
// // 			     _temperature[i] ,// the temperature field
// // 			     _actEnergy[i] ,
// // 			     _h->arrhA() ,
// // 			     R_g );
// // 	}
// // 	else{
// // 	  vc = creepvelocity(
// // 			     _thick,
// // 			     _tau[i] ,
// // 			     bkdT(int(i/_nd),i%_nd) ,// use background temp
// // 			     _actEnergy[i] ,
// // 			     _h->arrhA(), R_g );
// // 	}

// //         /* compute new loading rate*/
// //         slipdef_rate [i] = vp - vc;

// // 	/* adjust time step and record the max slip rate*/
// // 	//dt = min( fabs(umax / vc) , dtmax );
// // 	maxVc = max( abs(maxVc), vc );
// //       }
// //     dt = min( umax / maxVc, dtmax );

// //     /*     cout << "Max Slip = " << dt*maxVc << " m in "  */
// //     /* 	 << dt/second << " seconds (max = " << umax << " m)\n"; */

// //     return dt;
// //   }


// //   //============================================================
// //   /**
// //    * Get instantaneous creep deformation rate based on current
// //    * temperature and stress, return max slip rate (excluding
// //    * boundaries)
// //    **/
// //   double getSlipDefRate(
// // 			vector<double> &slipdef_rate // OUT: slip deficit rate
// // 			, const double &vp // plate velocity required to compute loading component
// // 			, const double &dtmax  // max allowed time step
// // 			, const double &maxslip  // max allowed slip in a time step
// // 			, double &t // current time used for temp computations
// // 			, bool isTemp // TRUE = consider temp changes, false uses background
// // 			 ) // INOUT current slip deficit rate, CHANGED here
// //   {

// //     vector<double> slipdef_rate_p(_nl*_nd, 0.0); // predicted

// //     /* Compute the slip deficit based on the current fault conditions */
// //     // computeStress();
// //     double dt_p = getInstSlipDefRate( slipdef_rate_p , vp, dtmax, maxslip,
// // 				      isTemp );
// //     // cout << "Changed time step from " << timestep / (year) << " yr to ";

// //     if( isTemp ){
// //       /* get heat rate for current slip deficit rate and stress*/
// //       vector<double> heat_rate(_nl*_nd, 0.0);
// //       getHeatRate( vp , slipdef_rate_p, heat_rate  );

// //       /* Predict stress after half a timestep */
// //       load( slipdef_rate_p , 0.5*dt_p );
// //       computeStress();

// //       /* Predict tempertature after half a timestep from previous slip*/
// //       computeTemp( t + 0.5*dt_p );

// //       /* Add on temperature due to slip in half time */
// //       for (int i = 0; i <_nl*_nd; ++i)/* add on temp from last half time step */
// // 	{
// // 	  /* get heat rate at new stress */
// // 	  double slip_rate = vp - slipdef_rate_p[i];
// // 	  double heatrate2 = slip_rate * _tau[i] * cell_area();
// //  	  _temperature[i] += integral_green_heat(
// // 						 heat_rate[i], // heat rate(t1)
// // 						 heatrate2, // heat rate (t2)
// // 						 t, // t1
// // 						 t + 0.5*dt_p , // t2
// // 						 t + 0.5*dt_p, // tnow
// // 						 _h->specHeat(), _h->density() ,
// // 						 _h->coolDist(), _h->diffusivity() , // constants for halfspace
// // 						 _cellLength , _cellDepth , _thick );
// // 	}
// //     } // end temp calc
// //     else{
// //       /* Predict stress after half a timestep */
// //       load( slipdef_rate_p , 0.5*dt_p );

// //       computeStress();

// //     }
// //     /* ReCompute the slip deficit based on the new fault conditions */
// //     double dt = getInstSlipDefRate( slipdef_rate , vp, dtmax, maxslip, isTemp );

// //     //cout << timestep / (year) << " year to ";
// //     //???timestep =  min( timestep , maxslip / maxVc ); // Compute time step to according to allowed slip
// //     //cout << timestep / (year) << " yr\n";

// //     /* Put fault back to previous condition */
// //     load( slipdef_rate_p , -0.5*dt_p );
// //     computeStress();

// //     return dt;
// //   }

// //   //============================================================
// //   /**
// //    * Get the heat rate on fault given current slip deficit rate and
// //    * current stress on fault
// //    **/
// //   void getHeatRate(
// // 		   const double &vp , // IN: plate motion
// // 		   vector<double> &slip_deficit_rate  , // IN: current slip deficit rate
// // 		   vector<double> &heat_rate // OUT: rate of heat generation at current state of fault
// // 		   )
// //   {
// //     if( heat_rate.size() != unsigned(_nl*_nd)  ) heat_rate.resize( _nl*_nd ); // check correct size of vector

// //     for( int i = 0 ; i < _nl*_nd ; i++ )
// //       {
// //         double slip_rate = vp - slip_deficit_rate[i];
// //         heat_rate[i] = slip_rate * _tau[i] * cell_area(); // TODO is this correct... Matthias original had _tau * slip_def_rate
// //       }

// //   }

// //   //============================================================
// //   /**
// //    *  Halve the time step until we have less than the desired number
// //    *  of critical cells
// //    **/
// //   bool adjustTimeStep(
// // 		      const int & maxNhypo , // maximum number of hypocenters
// // 		      vector<double> &slip_def_rate  , // current slip deficit rate
// // 		      double & timestep // time step to check
// // 		      )
// //   {
// //     bool is_adjusted = false;

// //     /** add on the slip */
// //     for( int i = 0 ; i < _nl*_nd ; i++ )
// //       {
// // 	double slipdef = slip_def_rate[i]*timestep;
// // 	load( int(i/_nd) , i%_nd , slipdef  );  // load the fault
// //       }

// //     /* Check we don't have too many critical cells, remove slip if we do */
// //     int tmpi , tmpj;
// //     computeStress();

// //     while( is_critical(tmpi,tmpj) > maxNhypo ){
// //       timestep *= 0.5; // halve the time step
// //       for( int i = 0 ; i < _nl*_nd ; i++ )
// // 	{
// // 	  add_slip( int(i/_nd) , i%_nd ,  slip_def_rate[i]*timestep );  // unload the fault
// // 	}
// //       computeStress(); // check again
// //       is_adjusted = true;
// //     }

// //     /* remove the other half */
// //     for( int i = 0 ; i < _nl*_nd ; i++ )
// //       {
// // 	add_slip( int(i/_nd) , i%_nd ,  slip_def_rate[i] * timestep );  // unload the fault
// //       }

// //     computeStress();
// //     return is_adjusted;
// //   }


// //   //============================================================
// //   /**
// //    * Compute temp change from frictional heat due to slip in a cell ,
// //    * return temperature change
// //    **/
// //   double frictionHeat(
// // 		      const int & i , const int & j, // coords
// // 		      const double &u // amount of slip
// // 		       ){

// //     /* Heat generated by friction */
// //     double heat = taud(i,j) * u * cell_area();

// //     /* new temp ... */
// //     double dT = green_heat( heat, 0, 0, // heat, time of slip, time now
// // 			    _h->specHeat(), _h->density() ,
// // 			    _h->coolDist(), _h->diffusivity() , // constants for halfspace
// // 			    _cellLength , _cellDepth , _thick ); // volume of cell
// //     return dT; // return temp change
// //   }


// //   //============================================================
// //   /**
// //    * Add record of earthquake slip to cells
// //    **/
// //   void add_eqk_record(
// // 		      vector <CellSlipEvent *> *cellslips , // pointer to vector of cell slip events
// // 		      vector <int> *strikeIdx ,
// // 		      vector <int> *depthIdx ) // pointer to vectors of slip locations
// //   {
// //     for( unsigned int ic = 0 ; ic < (*cellslips).size() ; ic++ ){
// //       int il = (*strikeIdx)[ic];
// //       int id = (*depthIdx)[ic];
// //       _slipHistory[il*_nd+id].add( (*cellslips)[ic] );
// //     }
// //     return;
// //   }

// //   //============================================================
// //   /**
// //    *  Add record of creep slip to cells
// //    **/
// //   void add_creep_record(
// // 			vector <CellCreepEvent *> *creepslips , // pointer to vector of cell slip events
// // 			const double & time ) // the time at end of time step
// //   {
// //     if( (*creepslips).size() != unsigned ( _nl*_nd ) ){
// //       cerr << "ERROR: creep not recorded correctly.\n";
// //     }
// //     else{
// //       for( int ic = 0 ; ic < _nl*_nd ; ic++ ){
// // 	_creepHistory[ic].add( (*creepslips)[ic] );
// //       }
// //     }

// //     return;
// //   }

// //   //============================================================

// //   //====================
// //   /**
// //      Predict time to failure for entire fault
// //   */
// //   double predictTimeToFailure( vector<double> &slip_def_rate ) // array of loading rates for entire fault
// //   {
// //     double min_time=1/0.0;

// //     for( int i = 0 ; i < _nl*_nd ; i++ )
// //       {
// //         if( slip_def_rate[i] > 0 ){ // don't consider cases where creep exceeds or is equal to plate loading
// //           double time = predictTimeToFailure( int(i/_nd) , i%_nd , slip_def_rate[i] );
// //           if( time < min_time ) min_time = time;
// //         }
// //       }
// //     return min_time;
// //   }

// //   //============================================================
// //   /**
// //    * Computing  percentage of the failure envelope on the fault
// //    **/
// //   double getStrength(
// // 		     const int &il, const int &id , // strike and depth index
// // 		     const double &prop_strength , // proportion of strength desired
// // 		     const double &vp ) // plate motion
// //   {
// //     double tau01 = prop_strength*taus(il,id);
// //     double tau02 = prop_strength*stress4creep( vp , _thick, temperature(il,id) ,
// // 					       _actEnergy[il*_nd+id] ,
// // 					       _h->arrhA() , R_g );
// //     return min(tau01,tau02);
// //   }

// //   //============================================================
// //   /**
// //    * Computing a background level of stress based on the minimum of
// //    * the arrest stress and the creep level of stress, similar to YBZ
// //    * 1991 code (see below)
// //    **/
// //   double getBkdStress(
// // 		     const int &il, const int &id , // strike and depth index
// // 		     const double &prop_strength , // proportion of strength desired
// // 		     const double &vp ) // plate motion
// //   {
// //     double tau01 = prop_strength*taua(il,id);
// //     double tau02 = prop_strength*stress4creep( vp , _thick, temperature(il,id) ,
// // 					       _actEnergy[il*_nd+id] ,
// // 					       _h->arrhA() , R_g );
// //     return min(tau01,tau02);
// //   }

// //   //============================================================
// //   /**
// //    * Computing a strength profile comparative to YBZ 1991 code
// //    */
// //   double getInitStress(
// // 		       const int &id , // depth index
// // 		       const double &dtaumx , // usually the cohesion, the maximum possible stress drop
// // 		       const double &vp ) // plate velocity
// //   {
// //     int il = int(0.5*_nl); // halfway along fault
// //     double tau01 = taus(il,id) - 0.5*dtaumx;
// //     double tau02 = 0.95*stress4creep( vp , _thick, temperature(il,id) ,
// // 				      _actEnergy[il*_nd+id] , _h->arrhA() ,
// // 				      R_g );
// //     return min(tau01,tau02);
// //   }

// //   //============================================================
// //   /**
// //    * Clean memory by deleting old creep events
// //    */
// //   void clean_mem(
// // 		 const double &time // the time now
// // 		  )
// //   {
// //     for( int i = 0 ; i < _nl*_nd; i++ )
// //       _creepHistory[i].clean( time - _creeptcut );

// //   }



// // };

