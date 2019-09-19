#include <mntMat3x3.h>

template<class T> 
Matrix3x3<T> dot(const Matrix3x3<T> &a, const Matrix3x3<T> &b) {
    Matrix3x3<T> res(0);
    for (size_t j = 0; j < 3; ++j) {
        for (size_t i = 0; i < 3; ++i) {
            for (size_t k = 0; k < 3; ++k) {
                res(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return res;
}


template<class T> 
Vector3<T> dot(const Matrix3x3<T> &a, const Vector3<T> &b) {
    Vector3<T> res(0);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            res[i] += a(i, j) * b[j];
        }
    }
    return res;
}

template<class T> 
Vector3<T> dot(const Vector3<T> &b, const Matrix3x3<T> &a) {
    Vector3<T> res(0);
    for (size_t j = 0; j < 3; ++j) {
        for (size_t i = 0; i < 3; ++i) {
            res[j] += b[i] * a(i, j);
        }
    }
    return res;
}

template<class T>
Matrix3x3<T> transpose(const Matrix3x3<T> &a) {
    Matrix3x3<T> res;
    for (size_t j = 0; j < 3; ++j) {
        for (size_t i = 0; i < 3; ++i) {
            res(j, i) = a(i, j);
        }
    }
}


