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


/** 2x2 matrix class.
  
  Matrices are column majored, ie
  elements in a column are adjacent in memory.  This
  implementation can be used with Lapack and is compatible with the
  Fortran array layout, but not with the "C/C++" array layout.
  
  */

template <class T>
class Matrix2x2
{

private:
  
  /// Stores the matrix elements.
  Vector4<T> v_;
  
public:


  /** Constructor with no arguments.  
    */
  Matrix2x2();                             
  
  /** Create matrix. Elelemnts are set.
    @param e value of each element
    */
  Matrix2x2(const T& e);

  /** Negation operator (elementwise).  */
  Matrix2x2<T> operator-() {
    for (size_t i = 0; i < 2*2; ++i)
      v_[i] *= -1;
  }

  /** Add "f" to every matrix
      element. */
  
  Matrix2x2<T> operator+=(const T f) {
    for (size_t i = 0; i < 2*2; ++i)
      v_[i] += f;    
  }

  /** Subtract "f" from each matrix
      element. */
  Matrix2x2<T> operator-=(const T f) {
    for (size_t i = 0; i < 2*2; ++i)
      v_[i] -= f;    
  }

  /** Multiply each matrix element by the value
      of "f". */
  Matrix2x2<T> operator*=(const T f) {
    for (size_t i = 0; i < 2*2; ++i)
      v_[i] *= f;    
  }

  /** Divide each matrix element by "f". */
  Matrix2x2<T> operator/=(const T f) {
    for (size_t i = 0; i < 2*2; ++i)
      v_[i] /= f;    
  }

  /** Add another matrix "b". */
  Matrix2x2<T> operator+=(const Matrix2x2<T> &b) {
    for (size_t i = 0; i < 2*2; ++i)
      v_[i] += b[i];    
  }

  /** Subtract another matrix "b" from itself. */
  Matrix2x2<T> operator-=(const Matrix2x2<T> &b) {
    for (size_t i = 0; i < 2*2; ++i)
      v_[i] -= b[i];    
  }

  /** Multiply each matrix element (of the matrix on
    the left of "*=") by the corresponding element of the "b" matrix*/
  Matrix2x2<T> operator*=(const Matrix2x2<T> &b) {
    for (size_t i = 0; i < 2*2; ++i)
      v_[i] *= b[i];    
  }

  /** Divide each matrix element (of the matrix on
    the left of "*=") by the corresponding element of the "b" matrix.*/
  Matrix2x2<T> operator/=(const Matrix2x2<T> &b) {
    for (size_t i = 0; i < 2*2; ++i)
      v_[i] /= b[i];    
  }

  /** Assignment operator. (default) */
  Matrix2x2<T> & operator=(const T&);

  /** Flat indexing operator.
    @param i flat index
    @return reference
    */
  inline T& operator[](size_t i)
  {
#ifndef NO_ASSERT
    assert(i < 2*2);
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
    assert(i < 2*2);
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
    assert(i < 2);
    assert(j < 2);
#endif 
    return v_[j * 2 + i];
  }

  /** Indexing operator. 
    @param i the row number
    @param j the column number
    @return reference to an element pointed by i and j
   */
  inline const T& operator()(size_t i, size_t j) const
  {
#ifndef NO_ASSERT
    assert(i < 2);
    assert(j < 2);
#endif 
    // column major
    return v_(j * 2 + i);
  }

};                                                              

typedef Matrix2x2<double> Mat2x2;

/** Matrix-matrix multiplication. Not to be confused with a*b, the elementwise
 multiplication. 
@see operator*
@param a a matrix
@param b another matrix of shape b.size(0)=a.size(1)
@return matrix = a.b
*/
template<class T> 
Matrix2x2<T> dot(const Matrix2x2<T> &a, const Matrix2x2<T> &b);

/** Matrix-vector multiplication. 
 @param a a matrix
 @param b a vector of length b.size()=a.size(1) 
 @return vector = a.b
*/
template<class T> 
Vector2<T> dot(const Matrix2x2<T> &a, const Vector2<T> &b);

/** Vector-matrix multiplication. 
 @param b a vector
 @param a a matrix of shape a.size(0) = b.size()
 @return vector = $b^T.a$
*/
template<class T> 
Vector2<T> dot(const Vector2<T> &b, const Matrix2x2<T> &a);

/** Transpose matrix
  @param a - input matrix
  @return matrix
*/
template<class T>
Matrix2x2<T> transpose(const Matrix2x2<T> &a);


#endif /* MNT_MAT2x2 */
