/* -*-C++-*- */
/* Created by A. Pletzer, Oct 23 1998, 2002
 * Changes by Irek Szczesniak, July 20 2001:
 *   - switched to the STL vector,
 *   - defined a new class Vector
 * $Id: MvVector.h 542 2013-10-18 02:22:16Z pletzer $
 */

#ifndef _VECTORF_H_
#define _VECTORF_H_    

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
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

// My headers
#include "MvFunctors.h"

/** The vector class

 The class to represent vectors.  Elements are adjacent in memory. */

template<class T>
class Vector : public std::vector<T> {  
                       
public:

  /*:::::::::::::::*/
  /* Constructors  */
  /*:::::::::::::::*/

  /** Constructor with no arguments (must subsequently invoke "alloc" to allocate memory) */
  Vector(); 
  
  /** Constructor: create vector of "n" elements.
    
    @param n vector size */
  Vector(size_t n); 

  /** Constructor: create vector of "n" elements "e".
    
    @param n vector size 
    @param e value of each element */
  Vector(size_t n, T e);

  /** Constructor: create vector from C array.
    
    @param vBeg start address 
    @param vEnd one beyond last address */
  Vector(const T* vBeg, const T* vEnd);

  /** Constructor: create vector from C array (C++11).
    
    @param array C array */
  Vector(std::initializer_list<T> array);

  /** Copy constructor: elements are copied into a new vector. 
   @param otherVec vector to be copied */
  Vector(const Vector<T>& otherVec);

  /** Assignment operator: set all elements to "f". 
   @param f scalar
   @return vector instance */
  Vector<T> &operator=(T f);

  /** Add the value "f" to every element of the vector. 
   @param f scalar 
   @return vector */
  Vector<T> &operator+=(T f);

  /** Subtract the value "f" to every element of the
      vector. 
  @param f scalar 
  @return vector */
  Vector<T> &operator-=(T f);

  /** Multiply every element by the value "f". 
  @param f scalar
  @return vector */
  Vector<T> &operator*=(T f);

  /** Divide every element by the value "f". 
  @param f scalar
  @return vector */
  Vector<T> &operator/=(T f);

  /** Vector addition. 
   @param w vector
   @return vector = original vector incremented by w */
  Vector<T> &operator+=(const Vector<T> &w);

  /** Vector subtraction. 
   @param w vector
   @return vector = original vector decremented by w  */
  Vector<T> &operator-=(const Vector<T> &w);

  /** Elementwise multiplication. Multiply every vector element (of the vector
    on the left of "*=") by its counterpart in the "b" vector.
  @param w vector
  @return vector */
  Vector<T> &operator*=(const Vector<T> &w);

  /** Elementwise division. 
    Divide every vector element (of the vector on
    the left of "*=") by its counterpart in the "b" vector.
  @param w vector
  @return vector  */
  Vector<T> &operator/=(const Vector<T> &w);

  /** Uniform grid generation. Elements are set so as to 
      span uniformly from xmin to xmax, in increasing order.
      @param xmin minimum value
      @param xmax maximum value
  */
  void space(T xmin, T xmax);

  /** Uniform grid generation with unit increment. Similar to space except that
      the increment is one. For instance x.range(0) has elements 0, 1, ... (up
      to size x.size()-1).
      @see space
  */
  void range(T imin=0);

  /** Fill the vector with random numbers between 0 and 1 (two subsequent calls 
      to random will generate different elements).
  */
  void random(void);

  

  /** Return the index of the element closest to the left of/at position of elem. 
      This assumes that the sequence is monotonically increasing. The 
      min/max return values are 0 and size()-1, respectively.
   @param elem reference abscissa
   @see ket
  */
  size_t bra(T elem);

  /** Return the indices of the elements that are closest to the left of elem. 
      This assumes that the sequence is monotonically increasing. 
   @param elem reference abscissa
  */
  Vector<size_t> bra(const Vector<T> &elem);

  
  /** Return the index of the element closest to the right of elem. 
      This assumes that the sequence is monotonically increasing. The
      min/max return values are 0 and size()-1, respectively.      
   @param elem reference abscissa
   @see bra
  */
  size_t ket(T elem);

  /** Return the indices of the elements that are closest to the right of elem. 
      This assumes that the sequence is monotonically increasing. 
   @param elem reference abscissa
  */
  Vector<size_t> ket(const Vector<T> &elem);

  /** Resizes the vector to the "n" size. 
   @param n new size
  */
  Vector<T> &alloc(size_t n);

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

  Vector<T> operator()(const Vector<size_t> &I) const;
};

typedef Vector<double> Vec;
typedef Vector<size_t> Vec_int;
typedef Vector< std::complex<double> > Vec_cmplx;

/**@name Vector Functions
  These are global functions operating on or generating vectors.  */

//@{

/** Addition. 
 @param v a vector
 @param w another vector
 @return vector = v + w 
*/
template<class T>
Vector<T> operator+(const Vector<T> &v, const Vector<T> &w);

/** Subtraction. 
 @param v a vector
 @param w another vector
 @return vector = v - w 
*/
template<class T>
Vector<T> operator-(const Vector<T> &v, const Vector<T> &w);

/** Elementwise multiplication. Not to be confused with the dot-product "dot".
 @see dot
 @param v a vector
 @param w another vector
 @return vector = v*w (!= dot(v, w))
*/
template<class T>
Vector<T> operator*(const Vector<T> &v, const Vector<T> &w);

/** Elementwise division. 
 @param v a vector
 @param w another vector
 @return vector = v/w */
template<class T>
Vector<T> operator/(const Vector<T> &v, const Vector<T> &w);

/** (Left) addition with a scalar. This is equivalent to creating a vector filled 
    with "f" and adding "a" to it.
    @param f a scalar
    @param a a vector
    @return vector = f + a
*/
template<class T>
Vector<T> operator+(T f, const Vector<T> &a);

/** (Right) addition with a scalar. This is equivalent to creating a vector filled 
    with "f" and adding "a" to it.
    @param a a vector
    @param f a scalar
    @return vector = a + f
*/
template<class T>
Vector<T> operator+(const Vector<T> &a, T f);

/** Subtraction from a scalar. This is equivalent to creating a vector filled with
    "f" and subtracting "w" from it.
    @param f a scalar
    @param w a vector
    @return vector = f-w
*/
template<class T>
Vector<T> operator-(T f, const Vector<T> &w);

/** Subtraction of a scalar. 
    @param w a vector
    @param f a scalar
    @return vector = w-f
*/
template<class T>
Vector<T> operator-(const Vector<T> &w, T f);

/** Negative. 
    @param v a vector
    @return vector = -w
*/
template<class T>
Vector<T> operator-(const Vector<T> &v);

/** (Left) multiplication of a vector by the scalar "f". 
 @param f s scalar
 @param a a vector
 @return vector = f*a
*/
template<class T>
Vector<T> operator*(T f, const Vector<T> &a);

/** (Right) multiplication of a vector by the scalar "f". 
 @param a a vector
 @param f s scalar
 @return vector = f*a
*/
template<class T>
Vector<T> operator*(const Vector<T> &a, T f);

/** Elementwise division. 
    @param f a scalar
    @param a a vector
    @return vector whose elements are f / element of "a"
*/
template<class T>
Vector<T> operator/(T f, const Vector<T> &a);

/** Scalar product. This is equivalent to sum(v*w). Not to be confused with the
 elementwise product v*w. 
@param v a vector
@param w another vector
@return vector = v.w
*/
template<class T>
T dot(const Vector<T> &v, const Vector<T> &w);

/** Apply function "sin" to each element. 
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> sin(const Vector<T> &v);

/** Apply function "cos" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> cos(const Vector<T> &v);

/** Apply function "tan" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> tan(const Vector<T> &v);

/** Apply function "asin" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> asin(const Vector<T> &v);

/** Apply function "acos" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> acos(const Vector<T> &v);

/** Apply function "atan" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> atan(const Vector<T> &v);

/** Apply function "exp" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> exp(const Vector<T> &v);

/** Apply function "log" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> log(const Vector<T> &v);

/** Apply function "sqrt" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> sqrt(const Vector<T> &v);

/** Apply function "abs" to each element.
    @param v a vector
    @return vector 
   */
template<class T>
Vector<T> abs(const Vector<T> &v);

/** Get the real part of a complex vector 
    @param v a vector
    @return real part of the vector
*/
Vec real(const Vec_cmplx &v);

/** Get the imaginary part of a complex vector 
    @param v a vector
    @return imaginary part of the vector
*/
Vec imag(const Vec_cmplx &v);

/** Get the conjugate of a complex vector 
    @param v a vector
    @return conjugate of the vector
*/
Vec_cmplx conjug(const Vec_cmplx &v);


/** Apply function "abs" to each element.
    @param v a vector
    @return real vector 
*/
Vec realabs(const Vec_cmplx& v);

/** Apply function "pow" to each element.
    @param v a vector
    @param exp the exponent
    @return vector 
 */
template<class T>
Vector<T> pow(const Vector<T> &v, T exp);

/** Apply function "pow" to each element.
    @param v a vector
    @param exp the exponent
    @return vector 
 */
template<class T>
Vector<T> pow(const Vector<T> &v, int exp);

/** Uniform grid generation. Create a vector of size n whose elements uniformly span 
    from xmin to xmax. 
    @param xmin 
    @param xmax
    @param n
 */
template<class T>
Vector<T> space(T xmin, T xmax, size_t n=2);

/** Uniform grid generation with integer increment. Create vector of size nsize whose
    elements are imin, imin+1, ... imin+size()-1.
    @see space
    @param imin
    @param nsize
    @return vector = [imin, imin+1, ... imin+nsize-1]
 */
template<class T>
Vector<T> range(T imin, size_t nsize=2);

/** Return the maximum value of v. 
    @param v vector 
    @return scalar = max(v)
 */
template<class T>
T max(const Vector<T> &v);

/** Take the maxium of two vectors.
    @param v a vector
    @param w another vector
    @return vector with elements $v_i > w_i ? v_i: w_i$.
 */
template<class T>
Vector<T> max(const Vector<T> &v, const Vector<T> &w);

/** Take the maximum of a vector and a scalar.
    @param v a vector
    @param f a scalar    
    @return vector with elements $v_i > f ? v_i: f$.
 */
template<class T>
Vector<T> max(const Vector<T> &v, T f);

/** Take the maximum of a vector and a scalar.
    @param f a scalar    
    @param v a vector
    @return vector with elements $v_i > f ? v_i: f$.
 */
template<class T>
Vector<T> max(T f, const Vector<T> &v);

/** Return the minimum element of v.
    @param v vector 
    @return scalar = min(v)
  */
template<class T>
T min(const Vector<T> &v);

/** Take the minimum of two vectors.
    @param v a vector
    @param w another vector
    @return vector with elements $v_i < w_i ? v_i: w_i$.
 */
template<class T>
Vector<T> min(const Vector<T> &v, const Vector<T> &w);

/** Take the minimum of a vector and a scalar.
    @param v a vector
    @param f a scalar    
    @return vector with elements $v_i < f ? v_i: f$.
 */
template<class T>
Vector<T> min(const Vector<T> &v, T f);

/** Take the minimum of a vector and a scalar.
    @param f a scalar    
    @param v a vector
    @return vector with elements $v_i < f ? v_i: f$.
 */
template<class T>
Vector<T> min(T f, const Vector<T> &v);

/** Sum all elements. 
    @see dot
    @param v a vector
    @return scalar = contraction of v.
 */
template<class T>
T sum(const Vector<T> &v);

/** Vector concatenation. 
    @param v a vector
    @param w another vector
    @return vector = concatenation of v and w
 */
template<class T>
Vector<T> cat( const Vector<T> &v, const Vector<T> &w);


/** Print out. 
 *
 * @param s stream
 * @param v Vector
 */
template <class T>
std::ostream& operator<<(std::ostream& s, const Vector<T>& v);

/** Set the real part of a complex vector
 * @param v complex vector to be modified
 * @param rV real part of the vector
 * @note the vector v must be pre-allocated to the correct size
 */
void setReal(Vec_cmplx &v, const Vec& rV);

/** Set the imaginary part of a complex vector
 * @param v complex vector to be modified
 * @param iV imaginary part of the vector
 * @note the vector v must be pre-allocated to the correct size
 */
void setImag(Vec_cmplx &v, const Vec& iV);

/** Set the real and the imaginary parts of a complex vector
 * @param v complex vector to be modified
 * @param rV real part of the vector
 * @param iV imaginary part of the vector
 * @note the vector v must be pre-allocated to the correct size
 */
void setRealImag(Vec_cmplx &v, const Vec& rV, const Vec& iV);

/** Create complex vector out of two real vectors
 * @param rV real part of the vector
 * @param iV imaginary part of the vector
 * @return complex vector
 */
Vec_cmplx cmplx(const Vec& rV, const Vec& iV);


//@}


#endif /* _VECTORF_H_ */
