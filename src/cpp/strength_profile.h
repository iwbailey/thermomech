#ifndef STRENGTH_PROFILE_H_
#define STRENGTH_PROFILE_H_

#include <vector>

//------------------------------------------------------------------------------
/**
    Class for defining the static (frictional) strength vs depth
*/
class StaticStrengthProfile{

public:
  /* Constructor */
  StaticStrengthProfile(
			const double &cohesion, // frictional cohesion
			const double &coeffFriction, // static frictional coefficient
			const double &dsigmaNdz, // normal stress gradient with depth
			const double &dsigmaPdz ):  // pore presure gradient
    _C(cohesion), _mu(coeffFriction), _dsigmadz( dsigmaNdz - dsigmaPdz){}

  /* Constructor for effective normal stress */
  StaticStrengthProfile(
			const double &cohesion, // frictional cohesion
			const double &coeffFriction, // static frictional coefficient
			const double &dsigmadz // effective normal stress gradient with depth
		       ):
    _C(cohesion), _mu(coeffFriction), _dsigmadz( dsigmadz){}

  /*Return the strength at specified depth */
  double operator() ( const double &depth ) const{
    return _C + _mu*_dsigmadz*depth;
  }

  /* Define on a fault grid, column ordered */
  std::vector<double> faultgrid_co( const int &nL, // dimensions along strike
                                    const int &nD, // dimensions down dip
                                    const double &cellHeight  // height of each slip cell
                                    ) const;


 private:
  /* Properties needed to compute the strength profile */
  double _C; // Cohesion
  double _mu;  // static coefficient of friction
  double _dsigmadz; // Effective normal stress gradient with depth

};


#endif /*STRENGTH_PROFILE_H_*/
