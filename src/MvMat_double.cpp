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
ColMat<double> eye(size_t n){
  ColMat<double> a(n, n, 0.0);
  for (size_t i = 0; i < n; ++i) 
      a(i,i) = 1.0;
  return a;
}


/* Functions */


/**********************************************************************

 load matrix from file
    
***********************************************************************/

ColMat<double> load(const std::string& cfile) { 

  // copy to temp1 

  std::string temp1 = "temp1"; 
  std::string s1 = "cp " + cfile + " " + temp1 + '\n'; 
  system(s1.c_str()); 

  // get rid of lines starting with % (comments) 

  std::string temp2 = "temp2"; 
  std::string s2 = "sed 's;^%.*$;;g' " + temp1 + " > " + temp2 + '\n'; 
  system(s2.c_str()); 
  
  const char *c_temp2 = temp2.c_str(); 
  ifstream file(c_temp2); 

  unsigned int nr, nc; 
  file >> nr; 
  file >> nc; 


  ColMat<double> a(nr, nc, 0.); 

  for(unsigned int i = 0; i < nr; ++i){ 
    for(unsigned int j = 0; j < nc; ++j){ 
      file >> a(i,j); 
    } 
  } 

  std::string s3 = "rm " + temp1 + ' ' + temp2 + '\n'; 
  system(s3.c_str()); 

  return a; 
} 

// real matrix dot complex vector
Vec_cmplx dot(const Mat &a, const Vec_cmplx &b)
{
  size_t i, j;
  Vec_cmplx c(a.size(0), std::complex<double>(0., 0.));

  size_t n = c.size();
  size_t n1 = a.size(1);
  for (i = 0; i < n; ++i)
    for (j = 0; j < n1; ++j)
      c(i) += a(i,j) * b(j);

  return c;
}

