#ifndef MNT_DOTS
#define MNT_DOTS

#include <mntMat2x2.h>
#include <mntMat2x3.h>
#include <mntMat3x2.h>
#include <mntMat3x3.h>
#include <mntVec2.h>
#include <mntVec3.h>

/**
 * Function that compute the dot product and the transpose
 */

template<class T> 
Matrix2x2<T> dot(const Matrix2x2<T> &a, const Matrix2x2<T> &b);

template<class T> 
Vector2<T> dot(const Matrix2x2<T> &a, const Vector2<T> &b);

template<class T> 
Vector2<T> dot(const Vector2<T> &b, const Matrix2x2<T> &a);


template<class T> 
Matrix2x2<T> dot(const Matrix2x3<T> &a, const Matrix3x2<T> &b);

template<class T> 
Vector2<T> dot(const Matrix2x3<T> &a, const Vector3<T> &b);

template<class T> 
Vector3<T> dot(const Vector2<T> &b, const Matrix2x3<T> &a);


template<class T> 
Matrix2x2<T> dot(const Matrix3x2<T> &a, const Matrix2x3<T> &b);

template<class T> 
Vector3<T> dot(const Matrix3x2<T> &a, const Vector2<T> &b);

template<class T> 
Vector2<T> dot(const Vector3<T> &b, const Matrix3x2<T> &a);


template<class T> 
Matrix3x3<T> dot(const Matrix3x3<T> &a, const Matrix3x3<T> &b);

template<class T> 
Vector3<T> dot(const Matrix3x3<T> &a, const Vector3<T> &b);

template<class T> 
Vector3<T> dot(const Vector3<T> &b, const Matrix3x3<T> &a);

////////////////////////////////////////////////////////////////

template<class T>
Matrix2x2<T> transpose(const Matrix2x2<T> &a);

template<class T>
Matrix3x2<T> transpose(const Matrix2x3<T> &a);

template<class T>
Matrix2x3<T> transpose(const Matrix3x2<T> &a);

template<class T>
Matrix3x3<T> transpose(const Matrix3x3<T> &a);

#endif /* MNT_DOTS */
