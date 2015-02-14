#include <iostream>
#include <vector>
#include "fileio.h"

int main( int argc, char *argv[])
{
  std::cerr << "Program: " << argv[0] << " No. Args: " << argc << std::endl;
  std::cerr << "jlanfkjjnevib\n";

  std::vector<double> tmp(128*32,0.0);
  txtfile2faultvec_dbl( tmp, "stressdrops_frac.txt", 128, 32 );


  std::cerr << tmp.size() << std::endl;

  for( unsigned int k=0; k<tmp.size(); k++ )
     std::cout << k << " " << tmp[k] << "\n";
}
