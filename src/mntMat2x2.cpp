#include <mntMat2x2.h>

template<class T> 
Matrix2x2<T> dot(const Matrix2x2<T> &a, const Matrix2x2<T> &b) {
    Matrix2x2<T> res(0);
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < 2; ++i) {
            for (size_t k = 0; k < 2; ++k) {
                res(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return res;
}


template<class T> 
Vector2<T> dot(const Matrix2x2<T> &a, const Vector2<T> &b) {
    Vector2<T> res(0);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            res[i] += a(i, j) * b[j];
        }
    }
    return res;
}

template<class T> 
Vector2<T> dot(const Vector2<T> &b, const Matrix2x2<T> &a) {
    Vector2<T> res(0);
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < 2; ++i) {
            res[j] += b[i] * a(i, j);
        }
    }
    return res;
}

template<class T>
Matrix2x2<T> transpose(const Matrix2x2<T> &a) {
    Matrix2x2<T> res;
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < 2; ++i) {
            res(j, i) = a(i, j);
        }
    }
}

// double
template class Matrix2x2<double>;
template Vector2<double> dot(const Matrix2x2<double>&, const Vector2<double>&);

