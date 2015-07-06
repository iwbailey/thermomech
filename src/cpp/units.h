#ifndef UNITS_H_
#define UNITS_H_
/**
   Namespace defining commonly used units for conversion purposes.
   Code uses std SI units, e.g., m, s, K.

   Example
   To output a distance variable from the code in km, you would use
   cout << distance/km << "km\n";

   To import a km value into the code, you would use
   cin >> distance_km;
   distance = distance_km*km;
 */
namespace Units
{

static const double pico = 1.0E-12; // pico
static const double nano = 1.0E-9; // nano
static const double mm = 1.0E-3; // milimeter
static const double cm = 1.0E-2; // centimeter
static const double m = 1.; // meter
static const double km = 1000; // kilometer
static const double kilo = 1.0E3; // kilo
static const double mega = 1.0E6; // mega
static const double giga = 1.0E9; // giga
static const double second = 1.; // second
static const double minute = 60.; // minute
static const double hour = 60*60; // hour
static const double day = 60*60*24; // day
static const double year = 60*60*24*365.25; // year
static const double g = 1.0E-3; // gramm
static const double kg = 1.; // Kilogram
static const double mole = 1; // mole
static const double K = 1; // degree Kelvin
static const double N = kg*m/(second*second); // Newton
static const double J = N*m; // Joule
static const double kJ = 1000*J; // Kilo joule
static const double erg = 1.0E-7*J; // erg
static const double Hz = 1/(second); // Hertz
static const double Pa = N/(m*m); // Pascal
static const double MPa = 1.0E6*Pa; // Mega Pscal
static const double GPa = 1.0E9*Pa; // Giga Pscal
static const double bar = 1.0E5*Pa; // bar
static const double RAD = 1.; // Radiant
static const double DEG= 3.1415 / 180; // Degree
static const double C = 1; // Coulomb
static const double eCharge = 1.602E-19 * C; // Elementary Charge
}

#endif /*UNITS_H_*/
