#ifndef MNT_DOTS
#define MNT_DOTS

#include <mntMatMxN.h>
#include <mntVecN.h>

/**
 * Functions that compute the dot product and the transpose
 */

template<int M, int N, int K, class T>
MatMxN<M, K, T> dot(const MatMxN<M, N, T> &a, const MatMxN<N, K, T> &b);

template<int M, int N, class T>
VecN<M, T> dot(const MatMxN<M, N, T> &a, const VecN<N, T> &b);

template<int M, int N, class T>
VecN<N, T> dot(const VecN<M, T> &b, const MatMxN<M, N, T> &a);

template<int M, class T>
T dot(const VecN<M, T> &a, const VecN<M, T> &b);

////////////////////////////////////////////////////////////////

template<int M, int N, class T>
MatMxN<N, M, T> transpose(const MatMxN<M, N, T> &a);

#endif /* MNT_DOTS */
