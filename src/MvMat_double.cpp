/**
 * list of functions applying to double matrices only 
 * $Id: MvMat_double.cpp 541 2013-10-17 19:55:14Z pletzer $
 **/

#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <complex>
#include <string>
#include "MvMatrix.h"

/* Identity matrix */
ColMat<double> eye(std::size_t n){
  ColMat<double> a(n, n, 0.0);
  for (std::size_t i = 0; i < n; ++i) 
      a(i,i) = 1.0;
  return a;
}


/* Functions */


/**********************************************************************

 load matrix from file
    
***********************************************************************/

ColMat<double> load(const std::string& cfile) {

  // strip lines starting with % (comments), in memory -- no temp files, no shelling
  // out to cp/sed/rm (which also meant cfile was being pasted, unescaped, into a shell
  // command)
  ifstream file(cfile.c_str());
  std::ostringstream stripped;
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] != '%') {
      stripped << line << '\n';
    }
  }

  std::istringstream in(stripped.str());

  unsigned int nr, nc;
  in >> nr;
  in >> nc;


  ColMat<double> a(nr, nc, 0.);

  for(unsigned int i = 0; i < nr; ++i){
    for(unsigned int j = 0; j < nc; ++j){
      in >> a(i,j);
    }
  }

  return a;
}

// real matrix dot complex vector
Vec_cmplx dot(const Mat &a, const Vec_cmplx &b)
{
  std::size_t i, j;
  Vec_cmplx c(a.size(0), std::complex<double>(0., 0.));

  std::size_t n = c.size();
  std::size_t n1 = a.size(1);
  for (i = 0; i < n; ++i)
    for (j = 0; j < n1; ++j)
      c(i) += a(i,j) * b(j);

  return c;
}

