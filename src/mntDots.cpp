#include "mntDots.h"

template<int M, int N, int K, class T> 
MatMxN<M, K, T> dot(const MatMxN<M, N, T> &a, const MatMxN<N, K, T> &b) {
    MatMxN<M, K, T> res( static_cast<T>(0) );
    for (int m = 0; m < M; ++m) {
        for (int k = 0; k < K; ++k) {
            for (int n = 0; n < N; ++n) {
                res(m, k) += a(m, n) * b(n, k);
            }
        }
    }
    return res;
}

template<int M, int N, class T> 
VecN<M, T> dot(const MatMxN<M, N, T> &a, const VecN<N, T> &b) {
    VecN<M, T> res( static_cast<T>(0) );
    for (int m = 0; m < M; ++m) {
        for (int n = 0; n < N; ++n) {
            res[m] += a(m, n) * b[n];
        }
    }
    return res;
}

template<int M, int N, class T> 
VecN<N, T> dot(const VecN<M, T> &b, const MatMxN<M, N, T> &a) {
    VecN<N, T> res( static_cast<T>(0) );
    for (int n = 0; n < N; ++n) {
        for (int  m = 0; m < M; ++m) {
            res[m] += b[n] * a(n, m);
        }
    }
    return res;
}


////////////////////////////////////////////////////////////////


template<int M, int N, class T>
MatMxN<N, M, T> transpose(const MatMxN<M, N, T> &a) {
    MatMxN<N, M, T> res;
    for (int m = 0; m < M; ++m) {
        for (int n = 0; n < N; ++n) {
            res(n, m) = a(m, n);
        }
    }
    return res;
}

////////////////////////////////////////////////////////////////////////
// template instantiations

// template MatMxN<2, 2, double> dot(const MatMxN<2, 2, double> &a, const MatMxN<2, 2, double> &b);
// template MatMxN<2, 3, double> dot(const MatMxN<2, 2, double> &a, const MatMxN<2, 3, double> &b);
// template MatMxN<3, 2, double> dot(const MatMxN<3, 2, double> &a, const MatMxN<2, 2, double> &b);
// template MatMxN<3, 3, double> dot(const MatMxN<3, 3, double> &a, const MatMxN<3, 3, double> &b);

// template VecN<2, double> dot(const MatMxN<2, 2, double> &a, const VecN<2, double> &b);
// template VecN<2, double> dot(const MatMxN<2, 3, double> &a, const VecN<3, double> &b);
// template VecN<3, double> dot(const MatMxN<3, 2, double> &a, const VecN<2, double> &b);
// template VecN<3, double> dot(const MatMxN<3, 3, double> &a, const VecN<3, double> &b);

// template VecN<2, double> dot(const VecN<2, double> &a, const MatMxN<2, 2, double> &b);
// template VecN<3, double> dot(const VecN<2, double> &a, const MatMxN<2, 3, double> &b);
// template VecN<2, double> dot(const VecN<3, double> &a, const MatMxN<3, 2, double> &b);
// template VecN<3, double> dot(const VecN<3, double> &a, const MatMxN<3, 3, double> &b);

// template MatMxN<2, 2, double> transpose(const MatMxN<2, 2, double> &a);
// template MatMxN<2, 3, double> transpose(const MatMxN<3, 2, double> &a);
// template MatMxN<3, 2, double> transpose(const MatMxN<2, 3, double> &a);
// template MatMxN<3, 3, double> transpose(const MatMxN<3, 3, double> &a);


