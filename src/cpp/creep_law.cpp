/**
    Source code for member functions of CreepLaw in creep_law.h
*/

#include "creep_law.h"

//------------------------------------------------------------------------------
/* Temperature dependent/stress independent part Arrhenius equation, i.e.
   A e^(-E/RT)
*/
double CreepLaw::temperatureTerm( const double &T )const
{
  return _ArrhAmplitude*exp( -_E/(PhysicalConstants::R_g * T ));
}
//------------------------------------------------------------------------------
/* Strain rate based on Arrhenius equation */
double CreepLaw::strainRate(
			    const double &tau, // stress
			    const double &T // temperature
			    ) const
{
  /* Stress dependency */
  double stressTerm = pow(tau, _n);

  /* Scale by amplitude to get strain rate */
  return  stressTerm * temperatureTerm(T);
}

//------------------------------------------------------------------------------
/* Stress for a given strain rate based on Arrhenius equation */
double CreepLaw::stress(
			const double &strainRate // desired creep rate
			, const double &T // temperature in K
			) const
{
  /* Scale the strain rate by the amplitude, combine temp dependence
     and apply the stress exponent scaling */
  return pow( strainRate/ temperatureTerm(T), 1.0/_n );
}

/* creep_law.cpp ends here */
