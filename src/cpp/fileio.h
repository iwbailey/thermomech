#ifndef FILEIO_H_
#define FILEIO_H_

#include <fstream>
#include <vector>
#include <string.h>
#include <iostream> // for runtime_error

/* Note templates functions must be defined where declared, can't be in separate
   .cpp file */

//------------------------------------------------------------------------------
/** Read a 2-d array of values as i, j, val columns of a text file into a
    vector */
template <typename T> bool txtfile2faultvec( std::vector<T> &outVec,
                                             const std::string &ifilename,
                                             const int &nL, const int &nD )
{
  /* Open the file as an ifstream */
  std::ifstream ifile( ifilename.c_str() );
  if( !ifile ){
    std::cerr << "Couldn't open the file " << ifilename <<  "\n";
    return false;
  }

  /* Define the output vector */
  outVec.resize(nL*nD);

  /* Define temp containers for the input */
  int i,j, nLines(0);
  while( ifile >> i )
    {
      ifile >> j;

      /* Read the value from the 3rd column, add it */
      ifile >> outVec[i*nD+j];
      nLines++;
    }

  /* Check total number of lines */
  if( nLines != nL*nD ){
    std::cerr << "Unexpected number of lines\n";
    return false;
  }

  /* Close the file */
  ifile.close();
  return true;
}

//------------------------------------------------------------------------------
/** Write a 2-d array of values in a vector into to a text file */
template <typename T> bool faultvec2txtfile( const std::string &ofilename,
                                             const std::vector<T> &faultvec,
                                             const int &nL, const int &nD )
{

  /* Open the file as an ofstream */
  std::ofstream ofile( ofilename.c_str() );

  //TODO: Check the file was successfully opened
  //TODO:  Check that the vector is the correct size

  /* Loop though the file in row order */
  for( unsigned int i=0; i<nL; i++ )
    for( unsigned int j=0; j<nD; j++ )
      {
        ofile << i << " " << j << " " << faultvec[i*nD+j] << std::endl;
      }

  /* Close the file */
  ofile.close();

  return true;
}

//------------------------------------------------------------------------------
#endif /* FILEIO_H_ */
