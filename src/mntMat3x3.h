#ifndef MNT_MAT3X3
#define MNT_MAT3X3

#ifndef NO_ASSERT
#include <cassert>
#endif

#include <math.h>
#include <stdlib.h>
#include <vector>
#include "mntVec3.h"
#include "mntVec9.h"


/** 3x3 matrix class.
  
  Matrices are column majored, ie
  elements in a column are adjacent in memory.  This
  implementation can be used with Lapack and is compatible with the
  Fortran array layout, but not with the "C/C++" array layout.
  
  */

template <class T>
class Matrix3x3
{

private:
  
  /// Stores the matrix elements.
  Vector9<T> v_;
  
public:


  /** Constructor with no arguments.  
    */
  Matrix3x3();                             
  
  /** Create matrix. Elelemnts are set.
    @param e value of each element
    */
  Matrix3x3(const T& e) {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] = e;
  }

  /** Negation operator (elementwise).  */
  Matrix3x3<T> operator-() {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] *= -1;
  }

  /** Add "f" to every matrix
      element. */
  
  Matrix3x3<T> operator+=(const T f) {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] += f;    
  }

  /** Subtract "f" from each matrix
      element. */
  Matrix3x3<T> operator-=(const T f) {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] -= f;    
  }

  /** Multiply each matrix element by the value
      of "f". */
  Matrix3x3<T> operator*=(const T f) {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] *= f;    
  }

  /** Divide each matrix element by "f". */
  Matrix3x3<T> operator/=(const T f) {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] /= f;    
  }

  /** Add another matrix "b". */
  Matrix3x3<T> operator+=(const Matrix3x3<T> &b) {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] += b[i];    
  }

  /** Subtract another matrix "b" from itself. */
  Matrix3x3<T> operator-=(const Matrix3x3<T> &b) {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] -= b[i];    
  }

  /** Multiply each matrix element (of the matrix on
    the left of "*=") by the corresponding element of the "b" matrix*/
  Matrix3x3<T> operator*=(const Matrix3x3<T> &b) {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] *= b[i];    
  }

  /** Divide each matrix element (of the matrix on
    the left of "*=") by the corresponding element of the "b" matrix.*/
  Matrix3x3<T> operator/=(const Matrix3x3<T> &b) {
    for (size_t i = 0; i < 3*3; ++i)
      v_[i] /= b[i];    
  }

  /** Assignment operator. (default) */
  Matrix3x3<T> & operator=(const T&);

  /** Flat indexing operator.
    @param i flat index
    @return reference
    */
  inline T& operator[](size_t i)
  {
#ifndef NO_ASSERT
    assert(i < 9);
#endif 
    return v_[i];
  }

  /** Flat indexing operator.
    @param i flat index
    @return reference
    */
  inline const T& operator[](size_t i) const
  {
#ifndef NO_ASSERT
    assert(i < 9);
#endif 
    return v_[i];
  }

  /** Indexing operator. 
    @param i the row number
    @param j the column number
    @return reference to an element pointed by i and j
   */
  inline T& operator()(size_t i, size_t j)
  {
#ifndef NO_ASSERT
    assert(i < 3);
    assert(j < 3);
#endif 
    return v_[j * 3 + i];
  }

  /** Indexing operator. 
    @param i the row number
    @param j the column number
    @return reference to an element pointed by i and j
   */
  inline const T& operator()(size_t i, size_t j) const
  {
#ifndef NO_ASSERT
    assert(i < 3);
    assert(j < 3);
#endif 
    // column major
    return v_(j * 3 + i);
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
