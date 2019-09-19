#include <mntMat3x3.h>

template<class T> 
Matrix3x3<T> dot(const Matrix3x3<T> &a, const Matrix3x3<T> &b) {
    Matrix3x3<T> res(0);
    for (size_t j = 0; j < MAT3X3_NDIMS; ++j) {
        for (size_t i = 0; i < MAT3X3_NDIMS; ++i) {
            for (size_t k = 0; k < MAT3X3_NDIMS; ++k) {
                res(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return res;
}


template<class T> 
Vector3<T> dot(const Matrix3x3<T> &a, const Vector3<T> &b) {
    Vector3<T> res(0);
    for (size_t i = 0; i < MAT3X3_NDIMS; ++i) {
        for (size_t j = 0; j < MAT3X3_NDIMS; ++j) {
            res[i] += a(i, j) * b[j];
        }
    }
    return res;
}

template<class T> 
Vector3<T> dot(const Vector3<T> &b, const Matrix3x3<T> &a) {
    Vector3<T> res(0);
    for (size_t j = 0; j < MAT3X3_NDIMS; ++j) {
        for (size_t i = 0; i < MAT3X3_NDIMS; ++i) {
            res[j] += b[i] * a(i, j);
        }
    }
    return res;
}

template<class T>
Matrix3x3<T> transpose(const Matrix3x3<T> &a) {
    Matrix3x3<T> res;
    for (size_t j = 0; j < MAT3X3_NDIMS; ++j) {
        for (size_t i = 0; i < MAT3X3_NDIMS; ++i) {
            res(j, i) = a(i, j);
        }
    }
    return res;
}

// double
template class Matrix3x3<double>;
template Vector3<double> dot(const Matrix3x3<double>&, const Vector3<double>&);
template Vector3<double> dot(const Vector3<double>&, const Matrix3x3<double>&);
template Matrix3x3<double> dot(const Matrix3x3<double>&, const Matrix3x3<double>&);
template Matrix3x3<double> transpose(const Matrix3x3<double>&);



