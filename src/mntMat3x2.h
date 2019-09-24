#ifndef MNT_MAT3x2
#define MNT_MAT3x2

#ifndef NO_ASSERT
#include <cassert>
#endif

#include <math.h>
#include <stdlib.h>
#include "mntVec6.h"

#define I_MAT3X2_NDIMS 3
#define J_MAT3X2_NDIMS 2

/** 3x2 matrix class.
  
  Matrices are column majored, ie
  elements in a column are adjacent in memory.  This
  implementation can be used with Lapack and is compatible with the
  Fortran array layout.
  
  */

template <class T>
class Matrix3x2 : public Vector6<T> {
  
public:

  /** Constructor with no arguments. */
  Matrix3x2() {};                           
  
  /** Create matrix. Elements are set.
    @param e value of each element
    */
  Matrix3x2(const T& e) {
    for (size_t i = 0; i < this->size(); ++i) {
      (*this)[i] = e;
    }
  }

  /** Indexing operator. 
    @param i the row number
    @param j the column number
    @return reference to an element pointed by i and j
   */
  inline T& operator()(size_t i, size_t j) {
#ifndef NO_ASSERT
    assert(i < I_MAT3X2_NDIMS);
    assert(j < J_MAT3X2_NDIMS);
#endif 
    return (*this)[j*I_MAT3X2_NDIMS + i];
  }

  /** Indexing operator. 
    @param i the row number
    @param j the column number
    @return reference to an element pointed by i and j
   */
  inline const T& operator()(size_t i, size_t j) const {
#ifndef NO_ASSERT
    assert(i < I_MAT3X2_NDIMS);
    assert(j < J_MAT3X2_NDIMS);
#endif 
    // column major
    return (*this)[j*I_MAT3X2_NDIMS + i];
  }

};                                                              

typedef Matrix3x2<double> Mat3x2;

#endif /* MNT_MAT3x2 */

