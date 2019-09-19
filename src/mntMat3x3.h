#ifndef MNT_MAT3x3
#define MNT_MAT3x3

#ifndef NO_ASSERT
#include <cassert>
#endif

#include <math.h>
#include <stdlib.h>
#include <vector>
#include "mntVec3.h"
#include "mntVec9.h"

#define MAT3X3_NDIMS 3


/** 3x3 matrix class.
  
  Matrices are column majored, ie
  elements in a column are adjacent in memory.  This
  implementation can be used with Lapack and is compatible with the
  Fortran array layout.
  */

template <class T>
class Matrix3x3 : public Vector9<T> {
  
public:

  /** Constructor with no arguments. */
  Matrix3x3() {};                           
  
  /** Create matrix. Elements are set.
    @param e value of each element
    */
  Matrix3x3(const T& e) {
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
    assert(i < MAT3X3_NDIMS);
    assert(j < MAT3X3_NDIMS);
#endif 
    return (*this)[j*MAT3X3_NDIMS + i];
  }

  /** Indexing operator. 
    @param i the row number
    @param j the column number
    @return reference to an element pointed by i and j
   */
  inline const T& operator()(size_t i, size_t j) const {
#ifndef NO_ASSERT
    assert(i < MAT3X3_NDIMS);
    assert(j < MAT3X3_NDIMS);
#endif 
    // column major
    return (*this)[j*MAT3X3_NDIMS + i];
  }

};                                                              

typedef Matrix3x3<double> Mat3x3;

/** Matrix-matrix multiplication. Not to be confused with a*b, the elementwise
 multiplication. 
  @see operator*
  @param a a matrix
  @param b another matrix of shape b.size(0)=a.size(1)
  @return matrix = a.b
*/
template<class T> 
Matrix3x3<T> dot(const Matrix3x3<T> &a, const Matrix3x3<T> &b);

/** Matrix-vector multiplication. 
  @param a a matrix
  @param b a vector of length b.size()=a.size(1) 
  @return vector = a.b
*/
template<class T> 
Vector3<T> dot(const Matrix3x3<T> &a, const Vector3<T> &b);

/** Vector-matrix multiplication. 
  @param b a vector
  @param a a matrix of shape a.size(0) = b.size()
  @return vector = $b^T.a$
*/
template<class T> 
Vector3<T> dot(const Vector3<T> &b, const Matrix3x3<T> &a);

/** Transpose matrix
  @param a - input matrix
  @return matrix
*/
template<class T>
Matrix3x3<T> transpose(const Matrix3x3<T> &a);


#endif /* MNT_MAT3x3 */

