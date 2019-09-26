#include "mntDots.h"
#include "mntMat2x2.h"
#include "mntMat2x3.h"
#include "mntMat3x2.h"
#include "mntMat3x3.h"

template<class T> 
Matrix2x2<T> dot(const Matrix2x2<T> &a, const Matrix2x2<T> &b) {
    Matrix2x2<T> res( static_cast<T>(0) );
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
    Vector2<T> res( static_cast<T>(0) );
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            res[i] += a(i, j) * b[j];
        }
    }
    return res;
}

template<class T> 
Vector2<T> dot(const Vector2<T> &b, const Matrix2x2<T> &a) {
    Vector2<T> res( static_cast<T>(0) );
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < 2; ++i) {
            res[j] += b[i] * a(i, j);
        }
    }
    return res;
}

////////////////////////////////////////////////////////////////

template<class T> 
Matrix2x2<T> dot(const Matrix2x3<T> &a, const Matrix3x2<T> &b) {
    Matrix2x2<T> res( static_cast<T>(0) );
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < 2; ++i) {
            for (size_t k = 0; k < 3; ++k) {
                res(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return res;
}


template<class T> 
Vector2<T> dot(const Matrix2x3<T> &a, const Vector3<T> &b) {
    Vector2<T> res( static_cast<T>(0) );
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            res[i] += a(i, j) * b[j];
        }
    }
    return res;
}

template<class T> 
Vector3<T> dot(const Vector2<T> &b, const Matrix2x3<T> &a) {
    Vector3<T> res( static_cast<T>(0) );
    for (size_t j = 0; j < 3; ++j) {
        for (size_t i = 0; i < 2; ++i) {
            res[j] += b[i] * a(i, j);
        }
    }
    return res;
}


////////////////////////////////////////////////////////////////


template<class T> 
Matrix3x3<T> dot(const Matrix3x2<T> &a, const Matrix2x3<T> &b) {
    Matrix3x3<T> res( static_cast<T>(0) );
    for (size_t j = 0; j < 3; ++j) {
        for (size_t i = 0; i < 3; ++i) {
            for (size_t k = 0; k < 2; ++k) {
                res(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return res;
}


template<class T> 
Vector3<T> dot(const Matrix3x2<T> &a, const Vector2<T> &b) {
    Vector3<T> res( static_cast<T>(0) );
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            res[i] += a(i, j) * b[j];
        }
    }
    return res;
}

template<class T> 
Vector2<T> dot(const Vector3<T> &b, const Matrix3x2<T> &a) {
    Vector2<T> res( static_cast<T>(0) );
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < 3; ++i) {
            res[j] += b[i] * a(i, j);
        }
    }
    return res;
}


////////////////////////////////////////////////////////////////

template<class T> 
Matrix3x3<T> dot(const Matrix3x3<T> &a, const Matrix3x3<T> &b) {
    Matrix3x3<T> res( static_cast<T>(0) );
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
    Vector3<T> res( static_cast<T>(0) );
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            res[i] += a(i, j) * b[j];
        }
    }
    return res;
}

template<class T> 
Vector3<T> dot(const Vector3<T> &b, const Matrix3x3<T> &a) {
    Vector3<T> res( static_cast<T>(0) );
    for (size_t j = 0; j < 3; ++j) {
        for (size_t i = 0; i < 3; ++i) {
            res[j] += b[i] * a(i, j);
        }
    }
    return res;
}


////////////////////////////////////////////////////////////////


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

template<class T>
Matrix3x2<T> transpose(const Matrix2x3<T> &a) {
    Matrix3x2<T> res;
    for (size_t j = 0; j < 3; ++j) {
        for (size_t i = 0; i < 2; ++i) {
            res(j, i) = a(i, j);
        }
    }
    return res;
}

template<class T>
Matrix2x3<T> transpose(const Matrix3x2<T> &a) {
    Matrix2x3<T> res;
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < 3; ++i) {
            res(j, i) = a(i, j);
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
    return res;
}

////////////////////////////////////////////////////////////////////////
// template instantiations


template Matrix2x2<double> dot(const Matrix2x2<double> &a, const Matrix2x2<double> &b);
template Vector2<double> dot(const Matrix2x2<double> &a, const Vector2<double> &b);
template Vector2<double> dot(const Vector2<double> &b, const Matrix2x2<double> &a);

template Matrix2x2<double> dot(const Matrix2x3<double> &a, const Matrix3x2<double> &b);
template Vector2<double> dot(const Matrix2x3<double> &a, const Vector3<double> &b);
template Vector3<double> dot(const Vector2<double> &b, const Matrix2x3<double> &a);

template Matrix3x3<double> dot(const Matrix3x2<double> &a, const Matrix2x3<double> &b);
template Vector3<double> dot(const Matrix3x2<double> &a, const Vector2<double> &b);
template Vector2<double> dot(const Vector3<double> &b, const Matrix3x2<double> &a);

template Matrix3x3<double> dot(const Matrix3x3<double> &a, const Matrix3x3<double> &b);
template Vector3<double> dot(const Matrix3x3<double> &a, const Vector3<double> &b);
template Vector3<double> dot(const Vector3<double> &b, const Matrix3x3<double> &a);


template Matrix2x2<double> transpose(const Matrix2x2<double> &a);
template Matrix3x2<double> transpose(const Matrix2x3<double> &a);
template Matrix2x3<double> transpose(const Matrix3x2<double> &a);
template Matrix3x3<double> transpose(const Matrix3x3<double> &a);

