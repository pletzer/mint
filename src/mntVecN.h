
#ifndef MNT_VECN
#define MNT_VECN 

// C headers
#include <cmath>
#include <cstdlib>
#include <cstring> // size_t
#include <ctime>
#ifndef NO_ASSERT
#include <cassert>
#endif

// STL
#include <array>
#include <numeric>
#include <algorithm>
#include <functional>

#include "MvFunctors.h"

/** Vector type allocated on the stack */

template<size_t N, class T>
class VecN : public std::array<T, N> {  
                       
public:

  /*:::::::::::::::*/
  /* Constructors  */
  /*:::::::::::::::*/

  /** Constructor with no arguments */
  VecN(); 
  
  /** Constructor: create vector of "e"'s.
    @param e value of each element */
  VecN(T e);

  /** Constructor: create vector from C pointer.
    @param ptr poointer */
  VecN(const T* ptr);

  /** Copy constructor: elements are copied into a new vector. 
   @param otherVec vector to be copied */
  VecN(const VecN<N, T>& otherVec);

  /** Assignment operator: set all elements to "f". 
   @param f scalar
   @return vector instance */
  VecN<N, T> &operator=(T f);

  /** Add the value "f" to every element of the vector. 
   @param f scalar 
   @return vector */
  VecN<N, T> &operator+=(T f);

  /** Subtract the value "f" to every element of the
      vector. 
  @param f scalar 
  @return vector */
  VecN<N, T> &operator-=(T f);

  /** Multiply every element by the value "f". 
  @param f scalar
  @return vector */
  VecN<N, T> &operator*=(T f);

  /** Divide every element by the value "f". 
  @param f scalar
  @return vector */
  VecN<N, T> &operator/=(T f);

  /** Vector addition. 
   @param w vector
   @return vector = original vector incremented by w */
  VecN<N, T> &operator+=(const VecN<N, T> &w);

  /** Vector subtraction. 
   @param w vector
   @return vector = original vector decremented by w  */
  VecN<N, T> &operator-=(const VecN<N, T> &w);

  /** Elementwise multiplication. Multiply every vector element (of the vector
    on the left of "*=") by its counterpart in the "b" vector.
  @param w vector
  @return vector */
  VecN<N, T> &operator*=(const VecN<N, T> &w);

  /** Elementwise division. 
    Divide every vector element (of the vector on
    the left of "*=") by its counterpart in the "b" vector.
  @param w vector
  @return vector  */
  VecN<N, T> &operator/=(const VecN<N, T> &w);

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

/**@name Vector Functions
  These are global functions operating on or generating vectors.  */

//@{

/** Addition. 
 @param v a vector
 @param w another vector
 @return vector = v + w 
*/
template<size_t N, class T>
VecN<N, T> operator+(const VecN<N, T> &v, const VecN<N, T> &w);

/** Subtraction. 
 @param v a vector
 @param w another vector
 @return vector = v - w 
*/
template<size_t N, class T>
VecN<N, T> operator-(const VecN<N, T> &v, const VecN<N, T> &w);

/** Elementwise multiplication. Not to be confused with the dot-product "dot".
 @see dot
 @param v a vector
 @param w another vector
 @return vector = v*w (!= dot(v, w))
*/
template<size_t N, class T>
VecN<N, T> operator*(const VecN<N, T> &v, const VecN<N, T> &w);

/** Elementwise division. 
 @param v a vector
 @param w another vector
 @return vector = v/w */
template<size_t N, class T>
VecN<N, T> operator/(const VecN<N, T> &v, const VecN<N, T> &w);

/** (Left) addition with a scalar. This is equivalent to creating a vector filled 
    with "f" and adding "a" to it.
    @param f a scalar
    @param a a vector
    @return vector = f + a
*/
template<size_t N, class T>
VecN<N, T> operator+(T f, const VecN<N, T> &a);

/** (Right) addition with a scalar. This is equivalent to creating a vector filled 
    with "f" and adding "a" to it.
    @param a a vector
    @param f a scalar
    @return vector = a + f
*/
template<size_t N, class T>
VecN<N, T> operator+(const VecN<N, T> &a, T f);

/** Subtraction from a scalar. This is equivalent to creating a vector filled with
    "f" and subtracting "w" from it.
    @param f a scalar
    @param w a vector
    @return vector = f-w
*/
template<size_t N, class T>
VecN<N, T> operator-(T f, const VecN<N, T> &w);

/** Subtraction of a scalar. 
    @param w a vector
    @param f a scalar
    @return vector = w-f
*/
template<size_t N, class T>
VecN<N, T> operator-(const VecN<N, T> &w, T f);

/** Negative. 
    @param v a vector
    @return vector = -w
*/
template<size_t N, class T>
VecN<N, T> operator-(const VecN<N, T> &v);

/** (Left) multiplication of a vector by the scalar "f". 
 @param f s scalar
 @param a a vector
 @return vector = f*a
*/
template<size_t N, class T>
VecN<N, T> operator*(T f, const VecN<N, T> &a);

/** (Right) multiplication of a vector by the scalar "f". 
 @param a a vector
 @param f s scalar
 @return vector = f*a
*/
template<size_t N, class T>
VecN<N, T> operator*(const VecN<N, T> &a, T f);

/** Elementwise division. 
    @param f a scalar
    @param a a vector
    @return vector whose elements are f / element of "a"
*/
template<size_t N, class T>
VecN<N, T> operator/(T f, const VecN<N, T> &a);

/** Scalar product. This is equivalent to sum(v*w). Not to be confused with the
 elementwise product v*w. 
@param v a vector
@param w another vector
@return vector = v.w
*/
template<size_t N, class T>
T dot(const VecN<N, T> &v, const VecN<N, T> &w);

/** Apply function "sin" to each element. 
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> sin(const VecN<N, T> &v);

/** Apply function "cos" to each element.
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> cos(const VecN<N, T> &v);

/** Apply function "tan" to each element.
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> tan(const VecN<N, T> &v);

/** Apply function "asin" to each element.
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> asin(const VecN<N, T> &v);

/** Apply function "acos" to each element.
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> acos(const VecN<N, T> &v);

/** Apply function "atan" to each element.
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> atan(const VecN<N, T> &v);

/** Apply function "exp" to each element.
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> exp(const VecN<N, T> &v);

/** Apply function "log" to each element.
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> log(const VecN<N, T> &v);

/** Apply function "sqrt" to each element.
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> sqrt(const VecN<N, T> &v);

/** Apply function "abs" to each element.
    @param v a vector
    @return vector 
   */
template<size_t N, class T>
VecN<N, T> abs(const VecN<N, T> &v);

/** Apply function "pow" to each element.
    @param v a vector
    @param exp the exponent
    @return vector 
 */
template<size_t N, class T>
VecN<N, T> pow(const VecN<N, T> &v, T exp);

/** Apply function "pow" to each element.
    @param v a vector
    @param exp the exponent
    @return vector 
 */
template<size_t N, class T>
VecN<N, T> pow(const VecN<N, T> &v, int exp);

/** Return the maximum value of v. 
    @param v vector 
    @return scalar = max(v)
 */
template<size_t N, class T>
T max(const VecN<N, T> &v);

/** Take the maxium of two vectors.
    @param v a vector
    @param w another vector
    @return vector with elements $v_i > w_i ? v_i: w_i$.
 */
template<size_t N, class T>
VecN<N, T> max(const VecN<N, T> &v, const VecN<N, T> &w);

/** Take the maximum of a vector and a scalar.
    @param v a vector
    @param f a scalar    
    @return vector with elements $v_i > f ? v_i: f$.
 */
template<size_t N, class T>
VecN<N, T> max(const VecN<N, T> &v, T f);

/** Take the maximum of a vector and a scalar.
    @param f a scalar    
    @param v a vector
    @return vector with elements $v_i > f ? v_i: f$.
 */
template<size_t N, class T>
VecN<N, T> max(T f, const VecN<N, T> &v);

/** Return the minimum element of v.
    @param v vector 
    @return scalar = min(v)
  */
template<size_t N, class T>
T min(const VecN<N, T> &v);

/** Take the minimum of two vectors.
    @param v a vector
    @param w another vector
    @return vector with elements $v_i < w_i ? v_i: w_i$.
 */
template<size_t N, class T>
VecN<N, T> min(const VecN<N, T> &v, const VecN<N, T> &w);

/** Take the minimum of a vector and a scalar.
    @param v a vector
    @param f a scalar    
    @return vector with elements $v_i < f ? v_i: f$.
 */
template<size_t N, class T>
VecN<N, T> min(const VecN<N, T> &v, T f);

/** Take the minimum of a vector and a scalar.
    @param f a scalar    
    @param v a vector
    @return vector with elements $v_i < f ? v_i: f$.
 */
template<size_t N, class T>
VecN<N, T> min(T f, const VecN<N, T> &v);

/** Sum all elements. 
    @see dot
    @param v a vector
    @return scalar = contraction of v.
 */
template<size_t N, class T>
T sum(const VecN<N, T> &v);

/** Print out. 
 *
 * @param s stream
 * @param v Vector
 */
template <size_t N, class T>
std::ostream& operator<<(std::ostream& s, const VecN<N, T>& v);

typedef VecN<2, double> Vec2;
typedef VecN<3, double> Vec3;
typedef VecN<4, double> Vec4;
typedef VecN<6, double> Vec6;
typedef VecN<9, double> Vec9;

//@}


#endif /* MNT_VECN */
