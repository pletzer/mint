/**
 * $Id: MvVector.cpp 548 2013-10-23 16:11:22Z pletzer $
 */

#include <complex>
#include "MvVector3.h"

/****************************************************************************/
// Definition of the Vector Class
/****************************************************************************/

template<class T>
Vector3<T>::Vector3()
{
}

template<class T>
Vector3<T>::Vector3(T e)
{
  for (size_t i = 0; i < this->size(); ++i)
  {
    (*this)[i] = e;
  }
}

template<class T>
Vector3<T>::Vector3(const Vector3<T> &w)
{
  (*this) = w;
}

template<class T>
void Vector3<T>::space(T xmin, T xmax)
{
  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = xmin + 
      (xmax-xmin)*static_cast<T>(i)/static_cast<T>(this->size()-1);
  }
}

template<class T> 
void Vector3<T>::random(void) 
{
  std::transform(this->begin(), this->end(), this->begin(), Random<T>());
}

template<class T>
void Vector3<T>::range(T imin) 
{
  for (size_t i = this->size(); i--;){
    (*this)[i] = imin + static_cast<T>(i);
  }
}

template<class T>
size_t Vector3<T>::bra(T elem)
{
  return static_cast<size_t>(
                 std::upper_bound(this->begin() + 1, 
                          this->end(), 
                          elem) - 
                 this->begin() - 1);
}

// upper_bound does not work for complex
template<>
size_t Vector3< std::complex<double> >::bra(std::complex<double> elem)
{
  size_t i = 0;
  while (std::abs(this->operator[](i)) < std::abs(elem)) {
    i++;
  }
  return i;
}

template<class T>
size_t Vector3<T>::ket(T elem){
  size_t n = this->size();
  size_t i = this->bra(elem) + 1;
  return ( (i<n)? i : n-1 );
}

template<class T>
Vector3<T> &Vector3<T>::operator=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), Setval<T>(f));
  return (*this);
}

template<class T>
Vector3<T> &Vector3<T>::operator+=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        std::bind2nd(std::plus<T>(), f));
  return (*this);
}

template<class T>
Vector3<T> &Vector3<T>::operator+=(const Vector3<T> &w) 
{
  std::transform(this->begin(), this->end(), 
             w.begin(), 
             this->begin(), 
             std::plus<T>() );
  return (*this);
}

template<class T>
Vector3<T> &Vector3<T>::operator-=(const Vector3<T> &w)
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::minus<T>() );
  return (*this);
}

template<class T>
Vector3<T> &Vector3<T>::operator-=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::minus<T>(), f));
  return (*this);
}

template<class T>
Vector3<T> &Vector3<T>::operator*=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::multiplies<T>(), f));
  return *this;  
}

template<class T>
Vector3<T> &Vector3<T>::operator*=(const Vector3<T> &w)
{  
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::multiplies<T>() );
  return (*this);  
}

template<class T>
Vector3<T> &Vector3<T>::operator/=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::divides<T>(), f));
  return (*this);
}

template<class T>
Vector3<T> &Vector3<T>::operator/=(const Vector3<T> &w) 
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::divides<T>() );
  return (*this);
}


template<class T>
Vector3<T> Vector3<T>::operator()(const Vector3<size_t> &I) const
{
  Vector3<T> y;

  for (size_t i = I.size(); i--;)
    y[i] = (*this)[ I[i] ];
  return y;
}

// Definition of extra functions supporting the Vector class

// find max element

template<class T>
T max(const Vector3<T> &v)
{
  return *std::max_element(v.begin(), v.end());
}

// find min element 

template<class T>
T min(const Vector3<T> &v) 
{
  return *std::min_element(v.begin(), v.end());
}

// the + operator

template<class T>
Vector3<T> operator+(const Vector3<T> &v, T f) 
{
  Vector3<T> c(v);
  return c += f;
}

template<class T>
Vector3<T> operator+(T f, const Vector3<T> &v) 
{
  Vector3<T> c(v);
  return c += f;
}

template<class T>
Vector3<T> operator+(const Vector3<T> &v, const Vector3<T> &w) 
{
  Vector3<T> c(v);
  return c += w;
}

template<class T>
Vector3<T> operator-(const Vector3<T> &v, T f)
{
  Vector3<T> c(v);
  return c -= f;
}

template<class T>
Vector3<T> operator-(T f, const Vector3<T> &v)
{
  Vector3<T> c = -v;
  return c += f;
}

template<class T>
Vector3<T> operator-(const Vector3<T> &v, const Vector3<T> &w)
{
  Vector3<T> c(v);
  return c -= w;
}

template<class T>
Vector3<T> operator-(const Vector3<T> &v)
{
  Vector3<T> c;
  std::transform(v.begin(), v.end(), c.begin(), std::negate<T>());
  return c;
}

template<class T>
T dot(const Vector3<T> &v, const Vector3<T> &w) 
{
  return std::inner_product(v.begin(), v.end(), w.begin(), static_cast<T>(0));
}

template<class T>
Vector3<T> operator*(const Vector3<T> &v, const Vector3<T> &w) 
{
  Vector3<T> c(v);
  return c*=w;
}

template<class T>
Vector3<T> operator*(T f, const Vector3<T> &v)
{
  Vector3<T> c(f);
  return c *= v;
}

template<class T>
Vector3<T> operator*(const Vector3<T> &v, T f)
{
  Vector3<T> c(f);
  return c *= v;
}

template<class T>
Vector3<T> operator/(const Vector3<T> &v, T f) 
{
  Vector3<T> c(v);
  return c/=f;
}

template<class T>
Vector3<T> operator/(T f, const Vector3<T> &v) 
{
  Vector3<T> c;
  std::transform(v.begin(), v.end(), c.begin(), bind1st(std::divides<T>(), f));
  return c;
}

// element by element division

template<class T>
Vector3<T> operator/(const Vector3<T> &v, const Vector3<T> &w) 
{
  Vector3<T> c(v);
  return c/=w;
}

// Math functions

template<class T>
Vector3<T> sin(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sin ) );
  return b;  
}  

template<class T>
Vector3<T> cos(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( cos ) );
  return b;  
}  

template<class T>
Vector3<T> tan(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( tan ) );
  return b;  
}  

template<class T>
Vector3<T> asin(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( asin ) );
  return b;  
}  

template<class T>
Vector3<T> acos(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( acos ) );
  return b;  
}  

template<class T>
Vector3<T> atan(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( atan ) );
  return b;  
}  

template<class T>
Vector3<T> exp(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( exp ) );
  return b;  
}  

template<class T>
Vector3<T> log(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( log ) );
  return b;  
}  

template<class T>
Vector3<T> sqrt(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sqrt ) );
  return b;  
}
  
template<class T>
Vector3<T> abs(const Vector3<T> &v) 
{ 
  Vector3<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), std::abs<T>);
  return b;  
}

template<>
Vector3< std::complex<double> > abs(const Vector3< std::complex<double> >& v) 
{
  Vector3< std::complex<double> > b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::abs(v[i]);
  return b;
}

Vec3 real(const Vec3_cmplx &v) {
  Vec3 b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::real(v[i]);
  return b;
}

Vec3 imag(const Vec3_cmplx &v) {
  Vec3 b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::imag(v[i]);
  return b;
}

Vec3_cmplx conjug(const Vec3_cmplx &v) {
  Vec3_cmplx b;
  for (size_t i = 0; i < v.size(); ++i) {
    b[i] = std::complex<double>(std::real(v[i]), -std::imag(v[i]));
  }
  return b;
}

Vec3 realabs(const Vec3_cmplx &v) 
{
  Vec3 b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::abs(v[i]);
  return b;
}

template<class T>
Vector3<T> pow(const Vector3<T> &v, T exp) 
{ 
  Vector3<T> b; 
  std::transform(v.begin(), v.end(), b.begin(), bind2nd(Power<T>(), exp));  
  return b;  
}  

template<class T>
Vector3<T> pow(const Vector3<T> &v, int exp) 
{ 
  // we've got to find a more efficient way to do this...
  return pow(v, static_cast<T>(exp));
}  

// grid generator
    
template<class T>
Vector3<T> space(T xmin, T xmax, size_t n)
{
  Vector3<T> a;
  for (size_t i = n; i--;)
    a[i] = xmin + (xmax - xmin)*T(i)/T(n - 1);
  return a;
}
    
template<class T>
Vector3<T> range(T imin, size_t n) 
{
  Vector3<T> a;
  for (size_t i = n; i--;)
    a[i] = imin + static_cast<T>(i);
  return a;
}

// max of 2 vectors

template<class T>
Vector3<T> max(const Vector3<T> &v, const Vector3<T> &w) 
{
  Vector3<T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) > w(i) ? v(i) : w(i);
  return res;
}

template<class T> 
Vector3<T> max(const Vector3<T> &v, T f) 
{
  Vector3<T> res;
  std::transform(v.begin(), v.end(), res.begin(), bind2nd(Max<T>(), f));
  return res;
}

template<class T>
Vector3<T> max(T f, const Vector3<T> &v) 
{
  return max(v, f);
}

// min of 2 arguments

template<class T> 
Vector3<T>  min(const Vector3<T> &v, const Vector3<T> &w) 
{
  Vector3<T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) < w(i) ? v(i) : w(i);
  return res;
}

template<class T> Vector3<T>  min(const Vector3<T> &v, T f) 
{
  Vector3<T> res;
  std::transform(v.begin(), v.end(), res.begin(), bind2nd(Min<T>(), f));
  return res;
}

template<class T>
Vector3<T> min(T f, const Vector3<T> &v) 
{
  return min(v, f);
}

/*
sum of elements
*/

template<class T>
T sum(const Vector3<T> &v) 
{
  return std::accumulate(v.begin(), v.end(), static_cast<T>(0));
}


// class methods to fill up elements of an existing Vector

template<class T>
void space(Vector3<T> &v, T xmin, T xmax) 
{
  size_t n = v.size();
  for (size_t i = n; i--;)
    v[i] = xmin + (xmax - xmin)*T(i)/T(n - 1);
}
    
template<class T>
void range(Vector3<T> &v, T imin) 
{
  for (size_t i = v.size(); i--;)
    v[i] = T(i) + imin;
}

// return index of min(v)

template<class T>
size_t index_min(Vector3<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = min_element(beg, end);
  return static_cast<size_t>(i - beg);
}

// return index of max(v)

template<class T>
size_t index_max(Vector3<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = max_element(beg, end);
  return static_cast<size_t>(i - beg);
}

void setReal(Vec3_cmplx &v, const Vec3& rV) {
  for (size_t i = 0; i < rV.size(); ++i) {
    v[i] = std::complex<double>(rV[i], v[i].imag());
  }
}

void setImag(Vec3_cmplx &v, const Vec3& iV) {
  for (size_t i = 0; i < iV.size(); ++i) {
    v[i] = std::complex<double>(v[i].real(), iV[i]);
  }
}

void setRealImag(Vec3_cmplx &v, const Vec3& rV, const Vec3& iV) {
  size_t n = rV.size() >= iV.size()? rV.size(): iV.size();
  for (size_t i = 0; i < n; ++i) {
    v[i] = std::complex<double>(rV[i], iV[i]);
  }
}

Vec3_cmplx cmplx(const Vec3& rV, const Vec3& iV) {
  Vec3_cmplx res;
  setRealImag(res, rV, iV);
  return res;
}

//*******************************************

template <class T>
std::ostream& operator<<(std::ostream& s, const Vector3<T>& v)
{
  for_each(v.begin(), v.end(), Print<T>(s));
  return s;
}

//
// template instantiations
//

// double
template class Vector3<double>;
template Vector3<double> space(double, double, size_t);
template double max(const Vector3<double>&);
template Vector3<double> max(const Vector3<double>&, double);
template Vector3<double> max(double, const Vector3<double>&);
template Vector3<double> max(const Vector3<double>&, const Vector3<double>&);
template double min(const Vector3<double>&);
template Vector3<double> min(const Vector3<double>&, double);
template Vector3<double> min(double, const Vector3<double>&);
template Vector3<double> min(const Vector3<double>&, const Vector3<double>&);
template double dot(const Vector3<double>&, const Vector3<double>&);
template double sum(const Vector3<double>&);
template Vector3<double> abs(const Vector3<double>&);
template Vector3<double> pow(const Vector3<double>&, double);
template Vector3<double> pow(const Vector3<double>&, int);
template Vector3<double> sqrt(const Vector3<double>&);
template Vector3<double> exp(const Vector3<double>&);
template Vector3<double> log(const Vector3<double>&);
template Vector3<double> sin(const Vector3<double>&);
template Vector3<double> cos(const Vector3<double>&);
template Vector3<double> tan(const Vector3<double>&);
template Vector3<double> asin(const Vector3<double>&);
template Vector3<double> acos(const Vector3<double>&);
template Vector3<double> atan(const Vector3<double>&);
template std::ostream& operator<<(std::ostream& s, const Vector3<double>&);
template Vector3<double> operator+(const Vector3<double>&, double);
template Vector3<double> operator+(double, const Vector3<double>&);
template Vector3<double> operator+(const Vector3<double>&, const Vector3<double>&);
template Vector3<double> operator-(const Vector3<double>&, double);
template Vector3<double> operator-(double, const Vector3<double>&);
template Vector3<double> operator-(const Vector3<double>&, const Vector3<double>&);
template Vector3<double> operator-(const Vector3<double>&);
template Vector3<double> operator*(const Vector3<double>&, double);
template Vector3<double> operator*(double, const Vector3<double>&);
template Vector3<double> operator*(const Vector3<double>&, const Vector3<double>&);
template Vector3<double> operator/(const Vector3<double>&, const Vector3<double>&);
template Vector3<double> operator/(const Vector3<double>&, double);
template Vector3<double> operator/(double, const Vector3<double>&);

// size_t
template class Vector3<size_t>;
template Vector3<size_t> range(size_t, size_t);
template Vector3<size_t> space(size_t, size_t, size_t);
template size_t max(const Vector3<size_t>&);
template Vector3<size_t> max(const Vector3<size_t>&, size_t);
template Vector3<size_t> max(size_t, const Vector3<size_t>&);
template Vector3<size_t> max(const Vector3<size_t>&, const Vector3<size_t>&);
template size_t min(const Vector3<size_t>&);
template Vector3<size_t> min(const Vector3<size_t>&, size_t);
template Vector3<size_t> min(size_t, const Vector3<size_t>&);
template Vector3<size_t> min(const Vector3<size_t>&, const Vector3<size_t>&);
template size_t dot(const Vector3<size_t>&, const Vector3<size_t>&);
template size_t sum(const Vector3<size_t>&);
template std::ostream& operator<<(std::ostream& s, const Vector3<size_t>&);
template Vector3<size_t> operator+(const Vector3<size_t>&, size_t);
template Vector3<size_t> operator+(size_t, const Vector3<size_t>&);
template Vector3<size_t> operator+(const Vector3<size_t>&, const Vector3<size_t>&);
template Vector3<size_t> operator-(const Vector3<size_t>&, size_t);
template Vector3<size_t> operator-(size_t, const Vector3<size_t>&);
template Vector3<size_t> operator-(const Vector3<size_t>&, const Vector3<size_t>&);
template Vector3<size_t> operator-(const Vector3<size_t>&);
template Vector3<size_t> operator*(const Vector3<size_t>&, size_t);
template Vector3<size_t> operator*(size_t, const Vector3<size_t>&);
template Vector3<size_t> operator*(const Vector3<size_t>&, const Vector3<size_t>&);

// int
template class Vector3<int>;
template Vector3<int> space(int, int, size_t);
template int max(const Vector3<int>&);
template Vector3<int> max(const Vector3<int>&, int);
template Vector3<int> max(int, const Vector3<int>&);
template Vector3<int> max(const Vector3<int>&, const Vector3<int>&);
template int min(const Vector3<int>&);
template Vector3<int> min(const Vector3<int>&, int);
template Vector3<int> min(int, const Vector3<int>&);
template Vector3<int> min(const Vector3<int>&, const Vector3<int>&);
template int dot(const Vector3<int>&, const Vector3<int>&);
template int sum(const Vector3<int>&);
template std::ostream& operator<<(std::ostream& s, const Vector3<int>&);
template Vector3<int> operator+(const Vector3<int>&, int);
template Vector3<int> operator+(int, const Vector3<int>&);
template Vector3<int> operator+(const Vector3<int>&, const Vector3<int>&);
template Vector3<int> operator-(const Vector3<int>&, int);
template Vector3<int> operator-(int, const Vector3<int>&);
template Vector3<int> operator-(const Vector3<int>&, const Vector3<int>&);
template Vector3<int> operator-(const Vector3<int>&);
template Vector3<int> operator*(const Vector3<int>&, int);
template Vector3<int> operator*(int, const Vector3<int>&);
template Vector3<int> operator*(const Vector3<int>&, const Vector3<int>&);


template class Vector3< std::complex<double> >;
template Vector3< std::complex<double> > space( std::complex<double> ,  std::complex<double> , size_t);
template std::complex<double>  dot(const Vector3< std::complex<double> >&, const Vector3< std::complex<double> >&);
template std::complex<double>  sum(const Vector3< std::complex<double> >&);
template Vector3< std::complex<double> > pow(const Vector3< std::complex<double> >&,  std::complex<double> );
template std::ostream& operator<<(std::ostream& s, const Vector3< std::complex<double> >&);
template Vector3< std::complex<double> > operator+(const Vector3< std::complex<double> >&,  std::complex<double> );
template Vector3< std::complex<double> > operator-(const Vector3< std::complex<double> >&,  std::complex<double> );
template Vector3< std::complex<double> > operator*(const Vector3< std::complex<double> >&,  std::complex<double> );
template Vector3< std::complex<double> > operator*( std::complex<double> , const Vector3< std::complex<double> >&);
template Vector3< std::complex<double> > operator+(const Vector3< std::complex<double> >&, const Vector3< std::complex<double> >&);
template Vector3< std::complex<double> > operator-(const Vector3< std::complex<double> >&, const Vector3< std::complex<double> >&);
template Vector3< std::complex<double> > operator*(const Vector3< std::complex<double> >&, const Vector3< std::complex<double> >&);
template Vector3< std::complex<double> > operator/(const Vector3< std::complex<double> >&, const Vector3< std::complex<double> >&);
template Vector3< std::complex<double> > operator-(const Vector3< std::complex<double> >&);


