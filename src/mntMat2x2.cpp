#include <mntMat2x2.h>

template<class T> 
Matrix2x2<T> dot(const Matrix2x2<T> &a, const Matrix2x2<T> &b) {
    Matrix2x2<T> res( static_cast<T>(0) );
    for (size_t j = 0; j < MAT2X2_NDIMS; ++j) {
        for (size_t i = 0; i < MAT2X2_NDIMS; ++i) {
            for (size_t k = 0; k < MAT2X2_NDIMS; ++k) {
                res(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return res;
}


template<class T> 
Vector2<T> dot(const Matrix2x2<T> &a, const Vector2<T> &b) {
    Vector2<T> res( static_cast<T>(0) );
    for (size_t i = 0; i < MAT2X2_NDIMS; ++i) {
        for (size_t j = 0; j < MAT2X2_NDIMS; ++j) {
            res[i] += a(i, j) * b[j];
        }
    }
    return res;
}

template<class T> 
Vector2<T> dot(const Vector2<T> &b, const Matrix2x2<T> &a) {
    Vector2<T> res( static_cast<T>(0) );
    for (size_t j = 0; j < MAT2X2_NDIMS; ++j) {
        for (size_t i = 0; i < MAT2X2_NDIMS; ++i) {
            res[j] += b[i] * a(i, j);
        }
    }
    return res;
}

template<class T>
Matrix2x2<T> transpose(const Matrix2x2<T> &a) {
    Matrix2x2<T> res;
    for (size_t j = 0; j < MAT2X2_NDIMS; ++j) {
        for (size_t i = 0; i < MAT2X2_NDIMS; ++i) {
            res(j, i) = a(i, j);
        }
    }
    return res;
}

// double
template class Matrix2x2<double>;
template Vector2<double> dot(const Matrix2x2<double>&, const Vector2<double>&);
template Vector2<double> dot(const Vector2<double>&, const Matrix2x2<double>&);
template Matrix2x2<double> dot(const Matrix2x2<double>&, const Matrix2x2<double>&);
template Matrix2x2<double> transpose(const Matrix2x2<double>&);



