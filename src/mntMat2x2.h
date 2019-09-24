#ifndef MNT_MAT2x2
#define MNT_MAT2x2

#ifndef NO_ASSERT
#include <cassert>
#endif

#include <math.h>
#include <stdlib.h>
#include <vector>
#include "mntVec2.h"
#include "mntVec4.h"

#define MAT2X2_NDIMS 2

/** 2x2 matrix class.
  
  Matrices are column majored, ie
  elements in a column are adjacent in memory.  This
  implementation can be used with Lapack and is compatible with the
  Fortran array layout.
  
  */

template <class T>
class Matrix2x2 : public Vector4<T> {
  
public:

  /** Constructor with no arguments. */
  Matrix2x2() {};                           
  
  /** Create matrix. Elements are set.
    @param e value of each element
    */
  Matrix2x2(const T& e) {
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
    assert(i < MAT2X2_NDIMS);
    assert(j < MAT2X2_NDIMS);
#endif 
    return (*this)[j*MAT2X2_NDIMS + i];
  }

  /** Indexing operator. 
    @param i the row number
    @param j the column number
    @return reference to an element pointed by i and j
   */
  inline const T& operator()(size_t i, size_t j) const {
#ifndef NO_ASSERT
    assert(i < MAT2X2_NDIMS);
    assert(j < MAT2X2_NDIMS);
#endif 
    // column major
    return (*this)[j*MAT2X2_NDIMS + i];
  }

};                                                              

typedef Matrix2x2<double> Mat2x2;

#endif /* MNT_MAT2x2 */

