#ifndef _CREEP_LAW_H
#define _CREEP_LAW_H
/**
   Class function defines a creep law for calculating aseismic slip
   rates on the fault
*/
#include <math.h>
#include "physical_constants.h"

//------------------------------------------------------------------------------
class CreepLaw{
  /**
     Creep law given by
     strainRate = A \tau^n exp{ -E/(R_g T) }

     where A = Arrhenius amplitude
     \tau = stress
     n = stress exponent
     E = activation energy
     R_g = Gas constant (defined in physical_constants.h)
     T = temperature
  */
 public:
  /**
     Constructor inputs are the Arrhenius eqn amplitude, the
     activation energy, and the stress exponent
  */
  CreepLaw( const double &ArrhAmplitude, const double &ActEnergy,
	    const double &stressExpon ) :
    _ArrhAmplitude(ArrhAmplitude), _E(ActEnergy), _n(stressExpon) {}

  /**
     Constructor inputs are the Arrhenius eqn amplitude, the
     activation energy. Uses default stress exponent of n=3
  */
  CreepLaw( const double &ArrhAmplitude, const double &ActEnergy ):
    _ArrhAmplitude(ArrhAmplitude), _E(ActEnergy), _n(3) {}

  /* Get strain rate for a given stress and temperature based on
     arrhenius eqn. Inputs are the stress and temperature */
  double strainRate( const double &tau, const double &T ) const;

  /* Get stress for a given strain rate based on arrhenius
     eqn. Inputs are the strain rate and temperature.*/
  double stress( const double &strainRate, const double &T ) const;

  /* Temperature dependent/stress independent part Arrhenius equation, i.e.
     A e^(-E/RT)
  */
  double temperatureTerm( const double &T ) const;

 private:
  /* Properties */
  double _ArrhAmplitude; // Arrhenius amplitude
  double _E; // Activation energy
  double _n; // Stress exponent
};

/* /\** */
/*  * Get the stress required for a certain aseismic slip rate */
/*  * Assume stress exponent is 3 */
/*  *\/ */

/* //###################################################################### */

/* /\** */
/*  * Get the activation energy required for a certain creep rate per Pa^3 */
/*  *\/ */
/* double acten4creep( */
/* 		   const double &vcPerPa3 // creep rate per Pa^3 */
/* 		   , const double &width // fault width in m */
/* 		   , const double &amp_arrh // arrhenius amplitude */
/* 		   , const double &R_g // gas constant */
/* 		   , const double &Temp // temperature, should use backgrnd */
/* 		   ) */
/* { */
/*   return - R_g * Temp * log( vcPerPa3 / (amp_arrh*width) ); */
/* } */
#endif
