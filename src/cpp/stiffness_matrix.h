#ifndef STIFFNESS_MATRIX
#define STIFFNESS_MATRIX
/**
   Class for stiffness matrix
*/

#include<math.h>
#include<vector>

class StiffnessMatrix{
 public:
  /* Constructor */
  StiffnessMatrix(
		  const double &dx, // cell length
		  const int &nx, // number cells length
		  const double &dz, // cell height
		  const int &nz, // number cells height
		  const double &rigidity // rgidity of halfspace
		  ): _nx(nx), _nz(nz), stiff_g(nx*nz*nz,0.0)
    {
      stiffness_init( dx, dz, 0.5*dz, rigidity );
    }

  /** Alternative constructor where we can specify upper slip cell depth */
  StiffnessMatrix(
		  const double &dx, // cell length
		  const int &nx, // number cells length
                  const double &z0, // depth of center of top slip cell
		  const double &dz, // cell height
		  const int &nz, // number cells height
		  const double &rigidity // rgidity of halfspace
		  ): _nx(nx), _nz(nz), stiff_g(nx*nz*nz,0.0)
    {
      stiffness_init( dx, dz, z0, rigidity );
    }

  /*Operator*/
  double operator() ( const int &k, const int &j1, const int &j2  ) const
  {
    return stiff_g[ j1*_nx*_nz + j2*_nx + k];
  }

 private:
  /* Properties */
  int _nx, _nz;
  std::vector<double> stiff_g;

  /* Allocate value to stiffness matrix */
  void set_stiff( const int &k, const int&j1, const int&j2, const double &val  ){
    stiff_g[ j1*_nx*_nz + j2*_nx + k] = val;
  }

  /* Initialize the stiffness matrix */
  void stiffness_init(
		      const double &dx,
		      const double &dz,
                      const double &z0,
		      const double &mu ); // rigidity

  /* Get the strain */
  double strn( const double &xw2, // Half width of a cell
	       const double &xd2, // Half height of a cell
	       const double &x3,  // depth of slipping cell center
	       const double &y1,  // distance along strike for stress calculation
	       const double &y3); // depth at which to calculate stress

};

#endif
