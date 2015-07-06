#ifndef FAULT_H_
#define FAULT_H_
/**
   Class definition for Fault object.  Holds grid of slip cells, their
   static and evolving properties.  Base class of ThermoFault
*/

#include <vector>
#include <algorithm>  // for max_element
#include <stdlib.h> // for abs
#include <stdexcept> // for runtime_error

#include "stiffness_matrix.h"
#include "creep_law.h"

//------------------------------------------------------------------------------
/**
   A fault is a grid of slip surfaces with no dip embedded in an
   elastic halfspace
*/
class Fault{

public:
  /** Constructor
  */
  /** Overloaded Constructor for a fault with no temperature decay **/
    Fault(
        const int &nL, // along-strike number of grid cells
        const int &nD, // down-dip number of grid cells
        const double &dx, // length of a single slip cell
        const double &dz, // height of a single slip cell
        const double &w, // fault width for strain rate calc
        const std::vector<double> &staticStrengths, // Static strength for each slip cell
        const std::vector<double> &dynamicStrengths, // dynamic strength for each slip cell
        const double &dynamicOvershootCoeff, // dynamic overshoot coefficient for fault
        const std::vector<double> &arrheniusAmpl, // amplitude in the arrhenius relation
        const std::vector<double> &stressExpon, // Stress exponent in the arrhenius relation
        const std::vector<double> &actEnergy, // activation energy
        const StiffnessMatrix &K, // Stiffness matrix for the fault
        const std::vector<double> &initStress // initial temperature for each slip cell
        ):
    _cellLength(dx), _cellHeight(dz), _nL( nL ), _nD( nD ), _faultWidth( w ),
    _Stiffness(K),
    _strengthStatic( staticStrengths ),
    _strengthDynamic( dynamicStrengths ),
    _dynamicOvershootCoeff( dynamicOvershootCoeff ),
    _arrhA( arrheniusAmpl ),
    _stressExpon( stressExpon ),
    _activEnergy( actEnergy ),
    _initStress(initStress),
    _slipDeficit(nL*nD, 0.0) {}

  //TODO: Check dimensions of the vectors in the constructor
public:

  /** Load the fault with constant plate motion for a fixed time period when
      creep velocity has already been calculated */
  void loadFault( const double &plateVelocity,
                  const std::vector<double> &creepVel, const double &dtime );

  /** Get the stress for a single slip cell at along strike coordinate i, depth
      j based on slip deficit accross the fault  */
  double getCellStress( const int &i, const int &j ) const;

  /** Get the stress for each slip cell on the fault */
  void getStress( std::vector<double> &stress ) const;

  /** Get the strain rate at a single slip cell */
  double getCellStrainRate( const int &k, const double &stress,
                            const double &temperature ) const;

  /** Get the creep velocity at a single slip cell  */
  double getCellCreepVelocity( const int &k, const double &stress,
                               const double &temperature ) const;

  /* ... overloaded version of above for i,j indexing */
  double getCellCreepVelocity( const int &i, const int &j, const double &stress,
                              const double &temperature )
  const { return getCellCreepVelocity( cellidx(i,j), stress, temperature ); }

  /** Get the creep velocity at all cells  */
  void getCreepVelocity( std::vector<double> &creepVel,
                         const std::vector<double> &stress,
                         const std::vector<double> &temperature ) const;

  /** Get the slip deficit of a single cell  */
  double getCellSlipDeficit( const int &i, const int &j )
  const { return _slipDeficit[cellidx(i,j)]; }

  /** Get the slip deficit of the entire fault */
  void getSlipDeficit( std::vector<double> &slipDef )
  const {
    slipDef = _slipDeficit;
  }

  /** Estimate time to failure for entire fault */
  double estimateTimeToFailure(
                               const std::vector<double> &stress,
                               const std::vector<double> &creepVel,
                               const double &plateVelocity ) const;

  /** Get the number of critical cells and coordinates of most critical cell
      given the current stress assuming static strength applies*/
  int nCriticalCells( unsigned int &iHypo, unsigned int&jHypo,
                      const std::vector<double> &stress ) const;

  /** Find the number of critical cells and their indices, incorporating dynamic
      strength */
  int findCriticalCells( std::vector<int> &criticalIdxs,
                                const std::vector<double> &stress,
                                const std::vector<bool> &isStatic ) const;


  /** Compute an earthquake given knowledge of currently critical hypocenter */
  void computeEarthquake( std::vector<double> &stress, // IN: current stress on fault, will be updated
                          double &time, // IN: current time, will be updated
                          const int&iHypo, const int &jHypo, // IN: coordinates of known hypocenter
                          const double &slipVelocity, // IN: velocity of all seismic slip
                          std::vector<double> &totalSlip ); // OUT: slip in each cell during the eqk

protected:
  /* Function used for consistent array indexing on the fault
     i = along strike index
     j = down-dip index
     out = index of vectors used to store fault properties
  */
  int cellidx( const int &i, const int &j ) const{ return i*_nD + j; }

  /** Get the arrest stress for cell k after failure*/
  double cellArrestStress( const int &k ) const{
    return _strengthStatic[k] -
      (_strengthStatic[k] - _strengthDynamic[k])*_dynamicOvershootCoeff;
  }

protected:
  /*Constant fault properties set at the begining*/
  int _nL;
  int _nD;
  double _cellLength, _cellHeight, _faultWidth;

  /*Fault dependencies*/
  StiffnessMatrix _Stiffness;
  double _dynamicOvershootCoeff;

  /* Constant cell properties, set at beginning*/
  std::vector <double>
    _initStress, // initial stress for each cell
    _strengthStatic, // static strength
    _strengthDynamic, // dynamic strength
    _arrhA, // arrhenius amplitude
    _stressExpon, // exponent to stress term in arrhenius relation
    _activEnergy; // Activation energy

  /* Dynamic cell properties*/
  std::vector <double> _slipDeficit;

};

//------------------------------------------------------------------------------
#endif /*FAULT_H_*/
 //  //============================================================
 // public:
 //  /**
 //     Initiation of stress; the static strength, background
 //   * temp and act energy must be initialized first
 //  */
 //  void initBkdStress(
 //                  const double & prop, // proportion of strength to set fault at
 //                  const double & vp ) // plate velocity
 //  {
 //    for( int i=0 ; i<_nl; i++ )
 //      for( int j=0 ; j<_nd ; j++ )
 //     _tau_bkd[i*_nd+j] = getBkdStress( i, j , prop, vp );

 //  }

 //  /* Alternate */
 //  void initBkdStress( const double & stress )// use a constant value
 //  {
 //    for( int i=0 ; i<_nl*_nd; i++ )
 //     _tau_bkd[i] = stress;
 //  }

 //  /* Alternate */
 //  void initBkdStress( vector<double> & stress )// use a constant value
 //  {
 //    for( int i=0 ; i<_nl*_nd; i++ )
 //     _tau_bkd[i] = stress[i];
 //  }

 //  /**
 //     Functions for access
 //  */
 //  int nl(){ return _nl; }  // get dimensions
 //  int nd(){ return _nd; }
 //  double cell_area(){ return _cellLength*_cellDepth; }
 //  double width(){ return _thick; }

 //  /**
 //      Get the background Temp at (i,j)
 //  */
 //  double bkdT( const int &i , const int & j ){
 //    double z = float((j+0.5)/_nd)*_depth;
 //    return (_h->Tsurface() +  _h->Tgrad() * z);
 //  }

 //  /**
 //     Get coordinates of center of cell
 //  */
 //  void getCoords( const int & i , const int & j , double & x , double & y ){
 //    x = float( (i+0.5)/_nl )*_length;
 //    y = float( (j+0.5)/_nd )*_depth;
 //  }
 //  double getDepth( const int & j ){ return float( (j+0.5)/_nd )*_depth; }

 //  /**
 //      Get background stress at (i,j)
 //  */
 //  double taub( const int &i , const int & j ){ return _tau_bkd[i*_nd+j]; }

 //  /**
 //     Get the static strength
 //  */
 //  double taus( const int &i , const int & j ){ return _taus[i*_nd+j]; }


 //  /**
 //      Get and set the arrest stress at (i,j)
 //  */
 //  double taua( const int &i , const int & j ){ return _taua[i*_nd+j]; }
 //  void set_taua( const int &i , const int & j , const double &taua ){
 //    _taua[i*_nd+j] = taua; }

 //  /**
 //      Get the dynamic strength at (i,j)
 //  */
 //  double taud( const int &i , const int & j ){
 //    return _taus[i*_nd+j] - ( _taus[i*_nd+j] - _taua[i*_nd+j] )/_doc ; }

 //  /**
 //     Get the activation energy at (i,j)
 //  */
 //  double actEnergy( const int &i , const int & j ){ return _actEnergy[i*_nd+j]; }

 //  void setActEnergy( vector<double> &actEn ){ initActEnergy( actEn ); }

 //  /**
 //     Get the current strength at (i,j)
 //  */
 //  double tauf( const int &i , const int & j ){ return (_tauf[i*_nd+j]);  }

 //  /**
 //      Get the stress at (i,j)
 //  */
 //  double stress( const int &i , const int & j ){ return _tau[i*_nd+j]; }

 //  /**
 //      Get the slip deficit at (i,j)
 //  */
 //  double slipdef( const int &i , const int & j ){ return _slipdef[i*_nd+j]; }

 //  /**
 //      Get the temp at (i,j)
 //  */
 //  double temperature( const int &i , const int & j ){ return _temperature[i*_nd+j]; }
 //  double T( const int &i , const int & j ){ return _temperature[i*_nd+j]; }

 //  /**
 //     Get the stress required for creep at vc
 //  */
 //  double creepstrength( const int &i, const int &j, const double &vc )
 //  {
 //    return stress4creep( vc , _thick, temperature(i,j) ,
 //                      _actEnergy[i*_nd+j] ,
 //                      _h->arrhA() , R_g );
 //  }

 //  /**
 //      Get the sitffness matrix value
 //  */
 //  double ff( const int &i , // difference in strike coords
 //          const int &j , // depth of cell where stress is computed (receiver)
 //          const int &k ){ // depth of slipping cell (source)
 //    return _stiffness[ k*_nl*_nd + j*_nl + i ];
 //  }

 //  /**
 //      Get the sitffness of a cell
 //  */
 //  double self_stiffness( const int &i , const int &j ){ // strike and depth indices
 //    return ff( 0 , j , j );
 //  }

 // /**
 //     Get Max stress
 // */
 //  double maxtau(){
 //    vector<double> :: iterator pos = max_element (_tau.begin(), _tau.end());
 //    return *pos;
 //  }

 // /**
 //     Get Max background stress
 // */
 //  double maxtaub(){
 //    vector<double> :: iterator pos = max_element (_tau_bkd.begin(), _tau_bkd.end());
 //    return *pos;
 //  }

 //  /**
 //     Get Max temperature
 //  */
 //  double maxTemp(){
 //    vector<double> :: iterator pos = max_element (
 //                                               _temperature.begin(),
 //                                               _temperature.end()
 //                                                );
 //    return *pos;
 //  }
 // /**
 //    Get Max background Temperature, used for setting up background stress
 // */
 //  double maxTbkd(){
 //    double T, maxT(0.0);
 //    for( int i = 0 ; i < _nl*_nd ; i++ ){
 //      T = bkdT(int(i/_nd),i%_nd);
 //      if( T > maxT ) maxT = T;
 //    }
 //    return maxT;
 //  }

 //  /**
 //      Get average temperature
 //  */
 //  double getAvgTemp(){
 //    double sumT = 0;
 //    for( int i = 0 ; i < _nl*_nd ; i++ ) sumT += _temperature[i];
 //    return sumT / (_nl*_nd);
 //  }

 //  /**
 //     Get average stress
 //  */
 //  double getAvgStress(){
 //    double sumS = 0;
 //    for( int i = 0 ; i < _nl*_nd ; i++ ) sumS += _tau[i];
 //    return sumS / (_nl*_nd);
 //  }

 //  /**
 //     Get average slip deficit
 //  */
 //  double getAvgSlipDef(){
 //    double sumSD = 0;
 //    for( int i = 0 ; i < _nl*_nd ; i++ ) sumSD += _slipdef[i];
 //    return sumSD / (_nl*_nd);
 //  }

 //  /**
 //     Get the seismic slip between the time limits t0 and t1
 //  */
 //  double getSeisSlip( const int &il , const int & id,
 //                   const double &t0, const double &t1 )
 //  {

 //    return _slipHistory[il*_nd+id].getSlip( t0, t1 );
 //  }

 //  /**
 //      Get the total seismic slip
 //  */
 //  double getTotalSeisSlip( const int &il , const int & id ) // coordinates
 //  {
 //    return _slipHistory[il*_nd+id].getSlip( 0.0 ); // get slip since time = 0.0
 //  }

 //  /**
 //     Get the total number of  seismic slip events
 //  */
 //  double getSeisSlipCount( const int &il , const int & id) // coordinates
 //  {
 //    return _slipHistory[il*_nd+id].getSlipCount();
 //  }

 //  /**
 //     Get the creep slip between the time limits t0 and t1
 //  */
 //  double getCreepSlip( const int &il , const int & id,
 //                   const double &t0, const double &t1 )
 //  {
 //    return _creepHistory[il*_nd+id].getSlip( t0, t1 );
 //  }

 //  /**
 //     Get the total aseismic slip
 //  */
 //  double getTotalCreepSlip( const int &il , const int & id ) // coordinates
 //  {
 //    return _creepHistory[il*_nd+id].getSlip( 0.0 ); // get slip since time = 0.0
 //  }


 //  /**
 //      A cell has failed so failure is now at dynamic levels
 //  */
 //  void become_dynamic( const int &i , const int & j ){ _tauf[i*_nd+j] = taud(i,j); }

 //  /**
 //   * Lock the fault
 //   */
 //  void lock(){ for( int i = 0 ; i < _nl*_nd ; i++ ) _tauf[i] = _taus[i]; }

 //  /**
 //   * Override the temperature - used in routines that check the change
 //   *  in strength due to temp
 //   */
 //  void setTemperature( const int &i , const int & j , const double & T ){
 //    _temperature[i*_nd+j] = T; }

 //  //============================================================
 //  /**
 //   * Return number of cells above failure threshold critical
 //   **/
 //  int is_critical(
 //               int &hi , int &hj ) // hypocenter coordinates
 //  {
 //    double maxtau = 0.0;
 //    int ncrit = 0;

 //    for( int il = 0 ; il < _nl ; il++ )
 //      for( int id = 0 ; id < _nd ; id++ ){ // loop over fault
 //     double overstress = _tau[il*_nd+id] - _tauf[il*_nd+id];
 //     if( overstress >= 0.0 ){ // check if critical
 //       ncrit++;
 //       if( overstress > maxtau ){  // check for hypocenter
 //         maxtau = overstress;
 //         hj = id; // record hypocenter
 //         hi = il;
 //       }
 //     }
 //      }

 //    return (ncrit);
 //  }

 //  //=======================
 //  /**
 //   * Return true if the locked fault is critical, return indices of
 //   * critical cells
 //   **/
 //  bool is_critical(
 //                vector<int> &crit_i ,
 //                vector<int> &crit_j ) //  coordinates of critical cells
 //  {
 //    crit_i.clear(); // empty containers
 //    crit_j.clear();
 //    bool isCritical = false;

 //    for( int il = 0 ; il < _nl ; il++ )
 //      for( int id = 0 ; id < _nd ; id++ ){ // loop over fault
 //     double overstress = _tau[il*_nd+id] - _tauf[il*_nd+id];
 //     if( overstress >= 0.0 ){ // check if critical
 //       isCritical = true;
 //       crit_i.push_back( il );
 //       crit_j.push_back( id );
 //     }
 //      }
 //    return isCritical;
 //  }

 //  //============================================================
 //  /**
 //   * Load the entire fault uniformly (except creepmask)
 //   */
 //  void load(
 //            const double &du )// amount of slip to add
 //  {
 //    for( int i = 0 ; i < _nl*_nd ; i++ )
 //      _slipdef[i] += _creepmask[i]*du;
 //  }

 //  /**
 //   * Load an individual cell
 //   **/
 //  void load( const int &i , const int & j , // coordinates
 //          double &delu )// INOUT: amount of slip to add
 //  {
 //    delu *= _creepmask[i*_nd+j]; // adjust for creep mask
 //    _slipdef[i*_nd+j] += delu;
 //  }

 //  /**
 //     load a fault given loading rate and time step
 //  */
 //  void load( vector<double> &loadrate , const double & dt )
 //  {
 //    for( int i = 0 ; i < _nl*_nd ; i++ )
 //      _slipdef[i] += _creepmask[i]*loadrate[i]*dt;
 //  }


 //  //============================================================
 //  /**
 //   *  Unload a cell
 //   **/
 //  void add_slip( const int &i , const int & j , // coordinates
 //              const double &delu ) // amount of slip deficit to remove
 //  { _slipdef[i*_nd+j] -= delu; }

 //  //============================================================
 //  /**
 //   * Overloaded copy of above Compute the stress from slip deficit and
 //   * return an array with stress values
 //   **/
 //  void computeStress( vector<double> &stress )
 //  {
 //    stress.clear();
 //    for (int i=0; i < _nl; ++i) // obs cell length index
 //      for (int j=0; j < _nd; ++j){ // obs cell depth index

 //     //_tau[i*_nd+j] = 0.0; // compute entirely from slip deficit
 //     _tau[i*_nd+j] = _tau_bkd[i*_nd+j]; // start with background level

 //     for (int l=0; l< _nd; ++l) // source cell depth index
 //       for (int k=0; k < _nl; ++k){ // source cell length index

 //         int ik = abs(int(i)- int(k)); //--- along strike difference
 //         _tau[i*_nd+j] += ff(ik,j,l) * ( -_slipdef[k*_nd+l] ); //  double checked order of indices!!
 //       } // end k and l
 //     stress.push_back(_tau[i*_nd+j]);

 //     if ( _tau[i*_nd+j] < 0.0 )
 //       cerr << "WARNING: -ve stress in cell [" << i << "," << j << "]\n";

 //      }// end i and j
 //    return;
 //  }


 //  /**
 //   * Compute the stress from slip deficit
 //   **/
 //  void computeStress()
 //  {
 //    vector<double> tmpstress;
 //    computeStress( tmpstress );
 //    return;
 //  }

 //  //============================================================

 //  /**
 //   * Get the temperature at (i,j) which is sum of background temp and
 //   * temp from heat generation on the fault Return average fault
 //   * temperature
 //   **/
 //  double computeTemp(
 //                  const double & time  // the time now
 //                  )
 //  {

 //    double oldT = getAvgTemp();

 //    for (int il=0; il < _nl; ++il)
 //      for (int id=0; id < _nd; ++id)
 //     {
 //       /* 2. add on the heat from the last q eqk events evaluated
 //          at the current time */
 //       double eqkT = _creepmask[il*_nd+id] *
 //         _slipHistory[il*_nd+id].getTempFromSlip(
 //                                                 time - _eqktcut ,
 //                                                 time ,
 //                                                 _cellLength ,
 //                                                 _cellDepth ,
 //                                                 _thick ,
 //                                                 taud(il,id) ,
 //                                                 _h );


 //       /* 3. add on the heat from the past creep */
 //       double creepT = _creepmask[il*_nd+id]*
 //         _creepHistory[il*_nd+id].getTempFromSlip(
 //                                                  time - _creeptcut ,
 //                                                  time ,
 //                                                  _cellLength ,
 //                                                  _cellDepth ,
 //                                                  _thick ,
 //                                                  _h );

 //       double Tij = bkdT(il,id) + creepT + eqkT;

 //       //if( Tij > bkdT(il,id) )
 //       //  cout << il << " " << id << " "
 //       //     << bkdT(il,id) << " " << creepT << " " << eqkT << endl;

 //       if( Tij < bkdT( il,id ) )
 //           cerr << "WARNING: Negative temp from slip at " << il << "," << id << endl;



 //       _temperature[il*_nd+id] = Tij;
 //     }

 //    //cerr << "Avg T change = " <<
 //    return getAvgTemp() - oldT;
 //  }

 //  //============================================================

 //  double computeActEnergy(
 //                       vector<double> &AE0  // Normal activation energy
 //                       , vector<double> &AE1  // Activation energy after eqk
 //                       , const double &t_h // healing time
 //                       , const double &t_now // time now
 //                        )
 //  {
 //    /**
 //       Compute activation energy based on time since last earthquake
 //     **/

 //    /* Loop through all cells */
 //    for( int il=0 ; il < _nl ; il++ )
 //      for( int id=0 ; id < _nd ; id++ )
 //     {
 //       /* Compute time since last eqk */
 //       double teq = _slipHistory[il*_nd+id].getLastSlipTime();

 //       if( teq < 0 || t_now - teq > t_h )
 //         /* Case where earthquake has no effect */
 //         _actEnergy[il*_nd+id ] = AE0[il*_nd+id];
 //       else if( teq == t_now  )
 //         _actEnergy[il*_nd+id ] = AE1[il*_nd+id];
 //       else{
 //         _actEnergy[il*_nd+id] = AE1[il*_nd+id] +
 //           (t_now - teq) * (AE0[il*_nd+id] - AE1[il*_nd+id])/t_h;
 //       }
 //     }
 //    return 0.0;
 //  }

 //  //============================================================
 //  /**
 //   * Distribute stress from change slip deficit
 //   **/
 //  void distributeStress( const int & ilc , const int &idc , const double delslip )
 //  {
 //    for (int il=0; il<_nl; ++il)
 //      for (int id=0; id<_nd; ++id)
 //     {
 //       int ik = abs(il - ilc);
 //       _tau[il*_nd+id] += ff(ik,id,idc) * delslip;
 //     }
 //  }

 //  //============================================================
 //  /**
 //   * Get creep vel based on current stress and temp
 //   **/
 //  void getInstSlipRate(
 //                    vector<double> &vc // OUT: column ordered vector
 //                    )
 //  {
 //    /* make sure container is right size */
 //    vc.resize( _nl*_nd );

 //    /* Compute the slip rate and time step */
 //    for (int i = 0; i <_nl*_nd; ++i)
 //      {
 //        /* compute creep velocity */
 //     vc[i] = creepvelocity( // see diffusion.h
 //                            _thick,
 //                            _tau[i] , // the stress field
 //                            _temperature[i] ,// the temperature field
 //                            _actEnergy[i] ,
 //                            _h->arrhA() ,
 //                            R_g );

 //      }

 //    return ;
 //  }


 //  //============================================================
 //  /**
 //   * Get instantaneous creep deformation rate based on
 //   * temperature and current stress, return time step
 //   **/
 //  double getInstSlipDefRate(
 //                         vector<double> &slipdef_rate // OUT: column ordered vector
 //                         , const double &vp  // plate velocity
 //                         , const double &dtmax // max time step
 //                         , const double &umax //max allowed creep
 //                         , bool isTemp // flag = true to use current temp, otherwise background
 //                         )
 //  {
 //    double dt; // timestep
 //    double maxVc( 0.0 );

 //    /* make sure container is right size */
 //    slipdef_rate.resize( _nl*_nd );

 //    /* Compute the slip rate and time step */
 //    for (int i = 0; i <_nl*_nd; ++i)
 //      {
 //        /* compute creep velocity */
 //     double vc;

 //     if( isTemp ){
 //       vc = creepvelocity( // see diffusion.h
 //                          _thick,
 //                          _tau[i] , // the stress field
 //                          _temperature[i] ,// the temperature field
 //                          _actEnergy[i] ,
 //                          _h->arrhA() ,
 //                          R_g );
 //     }
 //     else{
 //       vc = creepvelocity(
 //                          _thick,
 //                          _tau[i] ,
 //                          bkdT(int(i/_nd),i%_nd) ,// use background temp
 //                          _actEnergy[i] ,
 //                          _h->arrhA(), R_g );
 //     }

 //        /* compute new loading rate*/
 //        slipdef_rate [i] = vp - vc;

 //     /* adjust time step and record the max slip rate*/
 //     //dt = min( fabs(umax / vc) , dtmax );
 //     maxVc = max( abs(maxVc), vc );
 //      }
 //    dt = min( umax / maxVc, dtmax );

 //    /*     cout << "Max Slip = " << dt*maxVc << " m in "  */
 //    /*        << dt/second << " seconds (max = " << umax << " m)\n"; */

 //    return dt;
 //  }


 //  //============================================================
 //  /**
 //   * Get instantaneous creep deformation rate based on current
 //   * temperature and stress, return max slip rate (excluding
 //   * boundaries)
 //   **/
 //  double getSlipDefRate(
 //                     vector<double> &slipdef_rate // OUT: slip deficit rate
 //                     , const double &vp // plate velocity required to compute loading component
 //                     , const double &dtmax  // max allowed time step
 //                     , const double &maxslip  // max allowed slip in a time step
 //                     , double &t // current time used for temp computations
 //                     , bool isTemp // TRUE = consider temp changes, false uses background
 //                      ) // INOUT current slip deficit rate, CHANGED here
 //  {

 //    vector<double> slipdef_rate_p(_nl*_nd, 0.0); // predicted

 //    /* Compute the slip deficit based on the current fault conditions */
 //    // computeStress();
 //    double dt_p = getInstSlipDefRate( slipdef_rate_p , vp, dtmax, maxslip,
 //                                   isTemp );
 //    // cout << "Changed time step from " << timestep / (year) << " yr to ";

 //    if( isTemp ){
 //      /* get heat rate for current slip deficit rate and stress*/
 //      vector<double> heat_rate(_nl*_nd, 0.0);
 //      getHeatRate( vp , slipdef_rate_p, heat_rate  );

 //      /* Predict stress after half a timestep */
 //      load( slipdef_rate_p , 0.5*dt_p );
 //      computeStress();

 //      /* Predict tempertature after half a timestep from previous slip*/
 //      computeTemp( t + 0.5*dt_p );

 //      /* Add on temperature due to slip in half time */
 //      for (int i = 0; i <_nl*_nd; ++i)/* add on temp from last half time step */
 //     {
 //       /* get heat rate at new stress */
 //       double slip_rate = vp - slipdef_rate_p[i];
 //       double heatrate2 = slip_rate * _tau[i] * cell_area();
 //       _temperature[i] += integral_green_heat(
 //                                              heat_rate[i], // heat rate(t1)
 //                                              heatrate2, // heat rate (t2)
 //                                              t, // t1
 //                                              t + 0.5*dt_p , // t2
 //                                              t + 0.5*dt_p, // tnow
 //                                              _h->specHeat(), _h->density() ,
 //                                              _h->coolDist(), _h->diffusivity() , // constants for halfspace
 //                                              _cellLength , _cellDepth , _thick );
 //     }
 //    } // end temp calc
 //    else{
 //      /* Predict stress after half a timestep */
 //      load( slipdef_rate_p , 0.5*dt_p );

 //      computeStress();

 //    }
 //    /* ReCompute the slip deficit based on the new fault conditions */
 //    double dt = getInstSlipDefRate( slipdef_rate , vp, dtmax, maxslip, isTemp );

 //    //cout << timestep / (year) << " year to ";
 //    //???timestep =  min( timestep , maxslip / maxVc ); // Compute time step to according to allowed slip
 //    //cout << timestep / (year) << " yr\n";

 //    /* Put fault back to previous condition */
 //    load( slipdef_rate_p , -0.5*dt_p );
 //    computeStress();

 //    return dt;
 //  }

 //  //============================================================
 //  /**
 //   * Get the heat rate on fault given current slip deficit rate and
 //   * current stress on fault
 //   **/
 //  void getHeatRate(
 //                const double &vp , // IN: plate motion
 //                vector<double> &slip_deficit_rate  , // IN: current slip deficit rate
 //                vector<double> &heat_rate // OUT: rate of heat generation at current state of fault
 //                )
 //  {
 //    if( heat_rate.size() != unsigned(_nl*_nd)  ) heat_rate.resize( _nl*_nd ); // check correct size of vector

 //    for( int i = 0 ; i < _nl*_nd ; i++ )
 //      {
 //        double slip_rate = vp - slip_deficit_rate[i];
 //        heat_rate[i] = slip_rate * _tau[i] * cell_area(); // TODO is this correct... Matthias original had _tau * slip_def_rate
 //      }

 //  }

 //  //============================================================
 //  /**
 //   *  Halve the time step until we have less than the desired number
 //   *  of critical cells
 //   **/
 //  bool adjustTimeStep(
 //                   const int & maxNhypo , // maximum number of hypocenters
 //                   vector<double> &slip_def_rate  , // current slip deficit rate
 //                   double & timestep // time step to check
 //                   )
 //  {
 //    bool is_adjusted = false;

 //    /** add on the slip */
 //    for( int i = 0 ; i < _nl*_nd ; i++ )
 //      {
 //     double slipdef = slip_def_rate[i]*timestep;
 //     load( int(i/_nd) , i%_nd , slipdef  );  // load the fault
 //      }

 //    /* Check we don't have too many critical cells, remove slip if we do */
 //    int tmpi , tmpj;
 //    computeStress();

 //    while( is_critical(tmpi,tmpj) > maxNhypo ){
 //      timestep *= 0.5; // halve the time step
 //      for( int i = 0 ; i < _nl*_nd ; i++ )
 //     {
 //       add_slip( int(i/_nd) , i%_nd ,  slip_def_rate[i]*timestep );  // unload the fault
 //     }
 //      computeStress(); // check again
 //      is_adjusted = true;
 //    }

 //    /* remove the other half */
 //    for( int i = 0 ; i < _nl*_nd ; i++ )
 //      {
 //     add_slip( int(i/_nd) , i%_nd ,  slip_def_rate[i] * timestep );  // unload the fault
 //      }

 //    computeStress();
 //    return is_adjusted;
 //  }


 //  //============================================================
 //  /**
 //   * Compute temp change from frictional heat due to slip in a cell ,
 //   * return temperature change
 //   **/
 //  double frictionHeat(
 //                   const int & i , const int & j, // coords
 //                   const double &u // amount of slip
 //                    ){

 //    /* Heat generated by friction */
 //    double heat = taud(i,j) * u * cell_area();

 //    /* new temp ... */
 //    double dT = green_heat( heat, 0, 0, // heat, time of slip, time now
 //                         _h->specHeat(), _h->density() ,
 //                         _h->coolDist(), _h->diffusivity() , // constants for halfspace
 //                         _cellLength , _cellDepth , _thick ); // volume of cell
 //    return dT; // return temp change
 //  }


 //  //============================================================
 //  /**
 //   * Add record of earthquake slip to cells
 //   **/
 //  void add_eqk_record(
 //                   vector <CellSlipEvent *> *cellslips , // pointer to vector of cell slip events
 //                   vector <int> *strikeIdx ,
 //                   vector <int> *depthIdx ) // pointer to vectors of slip locations
 //  {
 //    for( unsigned int ic = 0 ; ic < (*cellslips).size() ; ic++ ){
 //      int il = (*strikeIdx)[ic];
 //      int id = (*depthIdx)[ic];
 //      _slipHistory[il*_nd+id].add( (*cellslips)[ic] );
 //    }
 //    return;
 //  }

 //  //============================================================
 //  /**
 //   *  Add record of creep slip to cells
 //   **/
 //  void add_creep_record(
 //                     vector <CellCreepEvent *> *creepslips , // pointer to vector of cell slip events
 //                     const double & time ) // the time at end of time step
 //  {
 //    if( (*creepslips).size() != unsigned ( _nl*_nd ) ){
 //      cerr << "ERROR: creep not recorded correctly.\n";
 //    }
 //    else{
 //      for( int ic = 0 ; ic < _nl*_nd ; ic++ ){
 //     _creepHistory[ic].add( (*creepslips)[ic] );
 //      }
 //    }

 //    return;
 //  }

 //  //============================================================
 //  /**
 //   * Compute time to failure given the current loading rate for a
 //   * single cell
 //   */
 //  double predictTimeToFailure( const int &il , const int & id , //coords
 //                            double slip_def_rate ) // loading rate for the cell
 //  {
 //    if( slip_def_rate < 0 )slip_def_rate = 0.0;
 //    double dtau = _taus[il*_nd+id] - _tau[il*_nd+id]; // stress before failure
 //    double du = -dtau/self_stiffness(il,id); // loading required
 //    return du/(_creepmask[il*_nd+id]*slip_def_rate);
 //  }

 //  //====================
 //  /**
 //     Predict time to failure for entire fault
 //  */
 //  double predictTimeToFailure( vector<double> &slip_def_rate ) // array of loading rates for entire fault
 //  {
 //    double min_time=1/0.0;

 //    for( int i = 0 ; i < _nl*_nd ; i++ )
 //      {
 //        if( slip_def_rate[i] > 0 ){ // don't consider cases where creep exceeds or is equal to plate loading
 //          double time = predictTimeToFailure( int(i/_nd) , i%_nd , slip_def_rate[i] );
 //          if( time < min_time ) min_time = time;
 //        }
 //      }
 //    return min_time;
 //  }

 //  //============================================================
 //  /**
 //   * Computing  percentage of the failure envelope on the fault
 //   **/
 //  double getStrength(
 //                  const int &il, const int &id , // strike and depth index
 //                  const double &prop_strength , // proportion of strength desired
 //                  const double &vp ) // plate motion
 //  {
 //    double tau01 = prop_strength*taus(il,id);
 //    double tau02 = prop_strength*stress4creep( vp , _thick, temperature(il,id) ,
 //                                            _actEnergy[il*_nd+id] ,
 //                                            _h->arrhA() , R_g );
 //    return min(tau01,tau02);
 //  }

 //  //============================================================
 //  /**
 //   * Computing a background level of stress based on the minimum of
 //   * the arrest stress and the creep level of stress, similar to YBZ
 //   * 1991 code (see below)
 //   **/
 //  double getBkdStress(
 //                  const int &il, const int &id , // strike and depth index
 //                  const double &prop_strength , // proportion of strength desired
 //                  const double &vp ) // plate motion
 //  {
 //    double tau01 = prop_strength*taua(il,id);
 //    double tau02 = prop_strength*stress4creep( vp , _thick, temperature(il,id) ,
 //                                            _actEnergy[il*_nd+id] ,
 //                                            _h->arrhA() , R_g );
 //    return min(tau01,tau02);
 //  }

 //  //============================================================
 //  /**
 //   * Computing a strength profile comparative to YBZ 1991 code
 //   */
 //  double getInitStress(
 //                    const int &id , // depth index
 //                    const double &dtaumx , // usually the cohesion, the maximum possible stress drop
 //                    const double &vp ) // plate velocity
 //  {
 //    int il = int(0.5*_nl); // halfway along fault
 //    double tau01 = taus(il,id) - 0.5*dtaumx;
 //    double tau02 = 0.95*stress4creep( vp , _thick, temperature(il,id) ,
 //                                   _actEnergy[il*_nd+id] , _h->arrhA() ,
 //                                   R_g );
 //    return min(tau01,tau02);
 //  }

 //  //============================================================
 //  /**
 //   * Clean memory by deleting old creep events
 //   */
 //  void clean_mem(
 //              const double &time // the time now
 //               )
 //  {
 //    for( int i = 0 ; i < _nl*_nd; i++ )
 //      _creepHistory[i].clean( time - _creeptcut );

 //  }


