#include <mntMatMxN.h>

template<size_t M, size_t K, size_t N, class T> 
MatMxN<M, N, T> dot(const MatMxN<M, K, T> &a, const MatMxN<K, N, T> &b) {
    MatMxN<M, N, T> res( static_cast<T>(0) );
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t k = 0; k < K; ++k) {
                res(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return res;
}

template<size_t M, size_t N, class T> 
VecN<M, T> dot(const MatMxN<M, N, T> &a, const VecN<N, T> &b) {
    VecN<M, T> res( static_cast<T>(0) );
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            res[i] += a(i, j) * b[j];
        }
    }
    return res;
}

template<size_t M, size_t N, class T> 
VecN<N, T> dot(const VecN<M, T> &b, const MatMxN<M, N, T> &a) {
    VecN<N, T> res( static_cast<T>(0) );
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < M; ++i) {
            res[j] += b[i] * a(i, j);
        }
    }
    return res;
}

template<size_t M, size_t N, class T>
MatMxN<N, M, T> transpose(const MatMxN<M, N, T> &a) {
    MatMxN<N, M, T> res;
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < M; ++i) {
            res(j, i) = a(i, j);
        }
    }
    return res;
}


// template instantiations
template class MatMxN<2, 2, double>;
template class MatMxN<2, 3, double>;
template class MatMxN<3, 2, double>;
template class MatMxN<3, 3, double>;

template VecN<2, double> dot(const MatMxN<2, 2, double> &a, const VecN<2, double> &b);
template VecN<2, double> dot(const MatMxN<2, 3, double> &a, const VecN<3, double> &b);
template VecN<3, double> dot(const MatMxN<3, 2, double> &a, const VecN<2, double> &b);
template VecN<3, double> dot(const MatMxN<3, 3, double> &a, const VecN<3, double> &b);

template VecN<2, double> dot(const VecN<2, double> &a, const MatMxN<2, 2, double> &b);
template VecN<3, double> dot(const VecN<2, double> &a, const MatMxN<2, 3, double> &b);
template VecN<2, double> dot(const VecN<3, double> &a, const MatMxN<3, 2, double> &b);
template VecN<3, double> dot(const VecN<3, double> &a, const MatMxN<3, 3, double> &b);

template MatMxN<2, 2, double> dot(const MatMxN<2, 2, double> &a, const MatMxN<2, 2, double> &b);
template MatMxN<2, 2, double> dot(const MatMxN<2, 3, double> &a, const MatMxN<3, 2, double> &b);
template MatMxN<2, 3, double> dot(const MatMxN<2, 2, double> &a, const MatMxN<2, 3, double> &b);
template MatMxN<2, 3, double> dot(const MatMxN<2, 3, double> &a, const MatMxN<3, 3, double> &b);
template MatMxN<3, 3, double> dot(const MatMxN<3, 3, double> &a, const MatMxN<3, 3, double> &b);

template MatMxN<2, 2, double> transpose(const MatMxN<2, 2, double> &a);
template MatMxN<3, 2, double> transpose(const MatMxN<2, 3, double> &a);
template MatMxN<2, 3, double> transpose(const MatMxN<3, 2, double> &a);
template MatMxN<3, 3, double> transpose(const MatMxN<3, 3, double> &a);



