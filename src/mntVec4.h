
#ifndef MNT_VEC4
#define MNT_VEC4

// C headers
#include <cmath>
#include <cstdlib>
#include <cstring> // size_t
#include <ctime>
#include <complex>
#ifndef NO_ASSERT
#include <cassert>
#endif

// STL
#include <array>
#include <numeric>
#include <algorithm>
#include <functional>

#include "MvFunctors.h"

/** The vector class

 The class to represent vectors.  Elements are adjacent in memory. */

template<class T>
class Vector4 : public std::array<T, 4> {  
                       
public:

  /*:::::::::::::::*/
  /* Constructors  */
  /*:::::::::::::::*/

  /** Constructor with no arguments */
  Vector4(); 
  
  /** Constructor: create vector of "n" elements "e".
    
    @param e value of each element */
  Vector4(T e);

  /** Copy constructor: elements are copied into a new vector. 
   @param otherVec vector to be copied */
  Vector4(const Vector4<T>& otherVec);

  /** Assignment operator: set all elements to "f". 
   @param f scalar
   @return vector instance */
  Vector4<T> &operator=(T f);

  /** Add the value "f" to every element of the vector. 
   @param f scalar 
   @return vector */
  Vector4<T> &operator+=(T f);

  /** Subtract the value "f" to every element of the
      vector. 
  @param f scalar 
  @return vector */
  Vector4<T> &operator-=(T f);

  /** Multiply every element by the value "f". 
  @param f scalar
  @return vector */
  Vector4<T> &operator*=(T f);

  /** Divide every element by the value "f". 
  @param f scalar
  @return vector */
  Vector4<T> &operator/=(T f);

  /** Vector addition. 
   @param w vector
   @return vector = original vector incremented by w */
  Vector4<T> &operator+=(const Vector4<T> &w);

  /** Vector subtraction. 
   @param w vector
   @return vector = original vector decremented by w  */
  Vector4<T> &operator-=(const Vector4<T> &w);

  /** Elementwise multiplication. Multiply every vector element (of the vector
    on the left of "*=") by its counterpart in the "b" vector.
  @param w vector
  @return vector */
  Vector4<T> &operator*=(const Vector4<T> &w);

  /** Elementwise division. 
    Divide every vector element (of the vector on
    the left of "*=") by its counterpart in the "b" vector.
  @param w vector
  @return vector  */
  Vector4<T> &operator/=(const Vector4<T> &w);

  /** Fill the vector with random numbers between 0 and 1 (two subsequent calls 
      to random will generate different elements).
  */
  void random(void);

  /** Sum of all elements 
   @return sum
  */
  T sum() const {
    T res = 0;
    for (size_t i = 0; i < this->size(); ++i) {
      res += (*this)(i);
    }
    return res;
  }

  /*::::::::::::::::::::::::::::::::*/                           
  /*  Index access operations */                           
  /*::::::::::::::::::::::::::::::::*/                           

  /** Return reference to "i"-th element.  
   @param i index 
   @return element */
  inline const T& operator()(size_t i) const
  {
#ifndef NO_ASSERT
    assert(i < this->size());
#endif
    return (*this)[i];
  }

  /** Returns reference to an element of index "i" (can be used for assignment). 
   @param i index 
   @return element */
  inline T& operator()(size_t i)
  {
#ifndef NO_ASSERT
    assert(i < this->size());
#endif
    return (*this)[i];
  }

};

typedef Vector4<double> Vec4;
typedef Vector4<size_t> Vec4_int;
typedef Vector4< std::complex<double> > Vec4_cmplx;

/**@name Vector Functions
  These are global functions operating on or generating vectors.  */

//@{

/** Addition. 
 @param v a vector
 @param w another vector
 @return vector = v + w 
*/
template<class T>
Vector4<T> operator+(const Vector4<T> &v, const Vector4<T> &w);

/** Subtraction. 
 @param v a vector
 @param w another vector
 @return vector = v - w 
*/
template<class T>
Vector4<T> operator-(const Vector4<T> &v, const Vector4<T> &w);

/** Elementwise multiplication. Not to be confused with the dot-product "dot".
 @see dot
 @param v a vector
 @param w another vector
 @return vector = v*w (!= dot(v, w))
*/
template<class T>
Vector4<T> operator*(const Vector4<T> &v, const Vector4<T> &w);

/** Elementwise division. 
 @param v a vector
 @param w another vector
 @return vector = v/w */
template<class T>
Vector4<T> operator/(const Vector4<T> &v, const Vector4<T> &w);

/** (Left) addition with a scalar. This is equivalent to creating a vector filled 
    with "f" and adding "a" to it.
    @param f a scalar
    @param a a vector
    @return vector = f + a
*/
template<class T>
Vector4<T> operator+(T f, const Vector4<T> &a);

/** (Right) addition with a scalar. This is equivalent to creating a vector filled 
    with "f" and adding "a" to it.
    @param a a vector
    @param f a scalar
    @return vector = a + f
*/
template<class T>
Vector4<T> operator+(const Vector4<T> &a, T f);

/** Subtraction from a scalar. This is equivalent to creating a vector filled with
    "f" and subtracting "w" from it.
    @param f a scalar
    @param w a vector
    @return vector = f-w
*/
template<class T>
Vector4<T> operator-(T f, const Vector4<T> &w);

/** Subtraction of a scalar. 
    @param w a vector
    @param f a scalar
    @return vector = w-f
*/
template<class T>
Vector4<T> operator-(const Vector4<T> &w, T f);

/** Negative. 
    @param v a vector
    @return vector = -w
*/
template<class T>
Vector4<T> operator-(const Vector4<T> &v);

/** (Left) multiplication of a vector by the scalar "f". 
 @param f s scalar
 @param a a vector
 @return vector = f*a
*/
template<class T>
Vector4<T> operator*(T f, const Vector4<T> &a);

/** (Right) multiplication of a vector by the scalar "f". 
 @param a a vector
 @param f s scalar
 @return vector = f*a
*/
template<class T>
Vector4<T> operator*(const Vector4<T> &a, T f);

/** Elementwise division. 
    @param f a scalar
    @param a a vector
    @return vector whose elements are f / element of "a"
*/
template<class T>
Vector4<T> operator/(T f, const Vector4<T> &a);

/** Scalar product. This is equivalent to sum(v*w). Not to be confused with the
 elementwise product v*w. 
@param v a vector
@param w another vector
@return vector = v.w
*/
template<class T>
T dot(const Vector4<T> &v, const Vector4<T> &w);

/** Apply function "sin" to each element. 
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> sin(const Vector4<T> &v);

/** Apply function "cos" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> cos(const Vector4<T> &v);

/** Apply function "tan" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> tan(const Vector4<T> &v);

/** Apply function "asin" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> asin(const Vector4<T> &v);

/** Apply function "acos" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> acos(const Vector4<T> &v);

/** Apply function "atan" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> atan(const Vector4<T> &v);

/** Apply function "exp" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> exp(const Vector4<T> &v);

/** Apply function "log" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> log(const Vector4<T> &v);

/** Apply function "sqrt" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> sqrt(const Vector4<T> &v);

/** Apply function "abs" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector4<T> abs(const Vector4<T> &v);

/** Get the real part of a complex vector 
    @param v a vector
    @return real part of the vector
*/
Vec4 real(const Vec4_cmplx &v);

/** Get the imaginary part of a complex vector 
    @param v a vector
    @return imaginary part of the vector
*/
Vec4 imag(const Vec4_cmplx &v);

/** Get the conjugate of a complex vector 
    @param v a vector
    @return conjugate of the vector
*/
Vec4_cmplx conjug(const Vec4_cmplx &v);


/** Apply function "pow" to each element.
    @param v a vector
    @param exp the exponent
    @return vector 
 */
template<class T>
Vector4<T> pow(const Vector4<T> &v, T exp);

/** Apply function "pow" to each element.
    @param v a vector
    @param exp the exponent
    @return vector 
 */
template<class T>
Vector4<T> pow(const Vector4<T> &v, int exp);

/** Return the maximum value of v. 
    @param v vector 
    @return scalar = max(v)
 */
template<class T>
T max(const Vector4<T> &v);

/** Take the maxium of two vectors.
    @param v a vector
    @param w another vector
    @return vector with elements $v_i > w_i ? v_i: w_i$.
 */
template<class T>
Vector4<T> max(const Vector4<T> &v, const Vector4<T> &w);

/** Take the maximum of a vector and a scalar.
    @param v a vector
    @param f a scalar    
    @return vector with elements $v_i > f ? v_i: f$.
 */
template<class T>
Vector4<T> max(const Vector4<T> &v, T f);

/** Take the maximum of a vector and a scalar.
    @param f a scalar    
    @param v a vector
    @return vector with elements $v_i > f ? v_i: f$.
 */
template<class T>
Vector4<T> max(T f, const Vector4<T> &v);

/** Return the minimum element of v.
    @param v vector 
    @return scalar = min(v)
  */
template<class T>
T min(const Vector4<T> &v);

/** Take the minimum of two vectors.
    @param v a vector
    @param w another vector
    @return vector with elements $v_i < w_i ? v_i: w_i$.
 */
template<class T>
Vector4<T> min(const Vector4<T> &v, const Vector4<T> &w);

/** Take the minimum of a vector and a scalar.
    @param v a vector
    @param f a scalar    
    @return vector with elements $v_i < f ? v_i: f$.
 */
template<class T>
Vector4<T> min(const Vector4<T> &v, T f);

/** Take the minimum of a vector and a scalar.
    @param f a scalar    
    @param v a vector
    @return vector with elements $v_i < f ? v_i: f$.
 */
template<class T>
Vector4<T> min(T f, const Vector4<T> &v);

/** Sum all elements. 
    @see dot
    @param v a vector
    @return scalar = contraction of v.
 */
template<class T>
T sum(const Vector4<T> &v);

/** Print out. 
 *
 * @param s stream
 * @param v Vector
 */
template <class T>
std::ostream& operator<<(std::ostream& s, const Vector4<T>& v);

/** Set the real part of a complex vector
 * @param v complex vector to be modified
 * @param rV real part of the vector
 */
void setReal(Vec4_cmplx &v, const Vec4& rV);

/** Set the imaginary part of a complex vector
 * @param v complex vector to be modified
 * @param iV imaginary part of the vector
 */
void setImag(Vec4_cmplx &v, const Vec4& iV);

/** Set the real and the imaginary parts of a complex vector
 * @param v complex vector to be modified
 * @param rV real part of the vector
 * @param iV imaginary part of the vector
 */
void setRealImag(Vec4_cmplx &v, const Vec4& rV, const Vec4& iV);

/** Create complex vector out of two real vectors
 * @param rV real part of the vector
 * @param iV imaginary part of the vector
 * @return complex vector
 */
Vec4_cmplx cmplx(const Vec4& rV, const Vec4& iV);


//@}


#endif /* MNT_VEC4 */
