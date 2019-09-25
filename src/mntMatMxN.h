#ifndef MNT_MATMxN
#define MNT_MATMxN

#ifndef NO_ASSERT
#include <cassert>
#endif

#include <cmath>
#include <cstdlib>
#include "mntVecN.h"


/** Column major matrix class */

template <size_t M, size_t N, class T>
class MatMxN : public VecN<N*M, T> {
  
public:

  /** Constructor with no arguments. */
  MatMxN() {};                           
  
  /** Create matrix. Elements are set.
    @param e value of each element
    */
  MatMxN(const T& e) {
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
    assert(i < M);
    assert(j < N);
#endif 
    return (*this)[i + j*M];
  }

  /** Indexing operator. 
    @param i the row number
    @param j the column number
    @return reference to an element pointed by i and j
   */
  inline const T& operator()(size_t i, size_t j) const {
#ifndef NO_ASSERT
    assert(i < M);
    assert(j < N);
#endif 
    // column major
    return (*this)[i + j*M];
  }

};

template<size_t M, size_t K, size_t N, class T> 
MatMxN<M, N, T> dot(const MatMxN<M, K, T> &a, const MatMxN<K, N, T> &b);

template<size_t M, size_t N, class T> 
VecN<M, T> dot(const MatMxN<M, N, T> &a, const VecN<N, T> &b);

template<size_t M, size_t N, class T> 
VecN<N, T> dot(const VecN<M, T> &b, const MatMxN<M, N, T> &a);

template<size_t M, size_t N, class T>
MatMxN<N, M, T> transpose(const MatMxN<M, N, T> &a);

typedef MatMxN<2, 2, double> Mat2x2;
typedef MatMxN<2, 3, double> Mat2x3;
typedef MatMxN<3, 2, double> Mat3x2;
typedef MatMxN<3, 3, double> Mat3x3;


#endif /* MNT_MATMxN */

