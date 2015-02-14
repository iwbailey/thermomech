#ifndef TEMPERATURE_PROFILE_H_
#define TEMPERATURE_PROFILE_H_

#include <vector>

//------------------------------------------------------------------------------
/**
    Class for the background temperature profile, defined so that we can quickly
    get the background temperature as a function of depth
 */
class BackgroundTemperatureProfile
{

 public:
  /*Constructor */
 BackgroundTemperatureProfile(
			      const double &surfaceTemperature, // Surface temperature
			      const double &temperatureGradient // Temperature gradient with depth
			      ):
   _T0(surfaceTemperature),
   _dTdz(temperatureGradient) {}

  /* Return the strength at specified depth */
  double operator() ( const double &depth ) const {
    return _T0 + _dTdz*depth;
  }

  /* Define on a fault grid, column ordered */
  std::vector<double> faultgrid_co( const int &nL, // number cells along strike
                                    const int &nD, // number cells down dip
                                    const double &cellHeight // height of each slip cell
                                    ) const;

 private:
  /* Properties needed to compute the background temperature profile */
  double _T0; // surface temperature in K
  double _dTdz;  // change in temperature with depth

};

#endif /*TEMPERATURE_PROFILE_H_*/
