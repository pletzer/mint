#ifndef MNT_MAT2x3
#define MNT_MAT2x3

#ifndef NO_ASSERT
#include <cassert>
#endif

#include <math.h>
#include <stdlib.h>
#include "mntMat2x2.h"
#include "mntMat3x2.h"
#include "mntMat3x3.h"
#include "mntVec2.h"
#include "mntVec3.h"
#include "mntVec6.h"

#define I_MAT2X3_NDIMS 2
#define J_MAT2X3_NDIMS 3

/** 2x3 matrix class.
  
  Matrices are column majored, ie
  elements in a column are adjacent in memory.  This
  implementation can be used with Lapack and is compatible with the
  Fortran array layout.
  
  */

template <class T>
class Matrix2x3 : public Vector6<T> {
  
public:

  /** Constructor with no arguments. */
  Matrix2x3() {};                           
  
  /** Create matrix. Elements are set.
    @param e value of each element
    */
  Matrix2x3(const T& e) {
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
    assert(i < I_MAT2X3_NDIMS);
    assert(j < J_MAT2X3_NDIMS);
#endif 
    return (*this)[j*I_MAT2X3_NDIMS + i];
  }

  /** Indexing operator. 
    @param i the row number
    @param j the column number
    @return reference to an element pointed by i and j
   */
  inline const T& operator()(size_t i, size_t j) const {
#ifndef NO_ASSERT
    assert(i < I_MAT2X3_NDIMS);
    assert(j < J_MAT2X3_NDIMS);
#endif 
    // column major
    return (*this)[j*I_MAT2X3_NDIMS + i];
  }

};                                                              

typedef Matrix2x3<double> Mat2x3;

#endif /* MNT_MAT2x3 */

