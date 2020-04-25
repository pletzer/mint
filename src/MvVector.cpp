/**
 * $Id: MvVector.cpp 548 2013-10-23 16:11:22Z pletzer $
 */

#include <complex>
#include "MvVector.h"

/****************************************************************************/
// Definition of the Vector Class
/****************************************************************************/

template<class T>
Vector<T>::Vector()
{
}

template<class T>
Vector<T>::Vector(size_t n) : std::vector<T>(n)
{
}

template<class T>
Vector<T>::Vector(size_t n, T e) : std::vector<T>(n, e)
{
}

template<class T>
Vector<T>::Vector(const T* vBeg, const T* vEnd) : std::vector<T>(vBeg, vEnd)
{
}

template<class T>
Vector<T>::Vector(std::initializer_list<T> array) : std::vector<T>{array}
{
}

template<class T>
Vector<T>::Vector(const Vector<T> &w)
{
  (*this) = w;
}

template<class T>
void Vector<T>::space(T xmin, T xmax)
{
  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = xmin + 
      (xmax-xmin)*static_cast<T>(i)/static_cast<T>(this->size()-1);
  }
}

template<class T> 
void Vector<T>::random(void) 
{
  std::transform(this->begin(), this->end(), this->begin(), Random<T>());
}

template<class T>
void Vector<T>::range(T imin) 
{
  for (size_t i = this->size(); i--;){
    (*this)[i] = imin + static_cast<T>(i);
  }
}

template<class T>
size_t Vector<T>::bra(T elem)
{
  return static_cast<size_t>(
                 std::upper_bound(this->begin() + 1, 
                          this->end(), 
                          elem) - 
                 this->begin() - 1);
}

// upper_bound does not work for complex
template<>
size_t Vector< std::complex<double> >::bra(std::complex<double> elem)
{
  size_t i = 0;
  while (std::abs(this->operator[](i)) < std::abs(elem)) {
    i++;
  }
  return i;
}

template<class T>
size_t Vector<T>::ket(T elem){
  size_t n = this->size();
  size_t i = this->bra(elem) + 1;
  return ( (i<n)? i : n-1 );
}

template<class T>
Vector<T> &Vector<T>::alloc(size_t n)
{
  this->resize(n);
  return (*this);
}

template<class T>
Vector<T> &Vector<T>::operator=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), Setval<T>(f));
  return (*this);
}

template<class T>
Vector<T> &Vector<T>::operator+=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        std::bind(std::plus<T>(), std::placeholders::_1, f));
  return (*this);
}

template<class T>
Vector<T> &Vector<T>::operator+=(const Vector<T> &w) 
{
#ifndef NO_ASSERT
  assert(this->size() == w.size());
#endif
  std::transform(this->begin(), this->end(), 
             w.begin(), 
             this->begin(), 
             std::plus<T>() );
  return (*this);
}

template<class T>
Vector<T> &Vector<T>::operator-=(const Vector<T> &w)
{
#ifndef NO_ASSERT
  assert(this->size() == w.size());
#endif
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::minus<T>() );
  return (*this);
}

template<class T>
Vector<T> &Vector<T>::operator-=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind(std::minus<T>(), std::placeholders::_1, f));
  return (*this);
}

template<class T>
Vector<T> &Vector<T>::operator*=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind(std::multiplies<T>(), std::placeholders::_1, f));
  return *this;  
}

template<class T>
Vector<T> &Vector<T>::operator*=(const Vector<T> &w)
{  
#ifndef NO_ASSERT
  assert(this->size() == w.size());
#endif
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::multiplies<T>() );
  return (*this);  
}

template<class T>
Vector<T> &Vector<T>::operator/=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        bind(std::divides<T>(), std::placeholders::_1, f));
  return (*this);
}

template<class T>
Vector<T> &Vector<T>::operator/=(const Vector<T> &w) 
{
#ifndef NO_ASSERT
  assert(this->size() == w.size());
#endif
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::divides<T>() );
  return (*this);
}


template<class T>
Vector<T> Vector<T>::operator()(const Vector<size_t> &I) const
{
  Vector<T> y(I.size());

  for (size_t i = I.size(); i--;)
    y[i] = (*this)[ I[i] ];
  return y;
}

// Definition of extra functions supporting the Vector class

// find max element

template<class T>
T max(const Vector<T> &v)
{
#ifndef NO_ASSERT
  assert(v.size()>0);
#endif
  return *std::max_element(v.begin(), v.end());
}

// find min element 

template<class T>
T min(const Vector<T> &v) 
{
#ifndef NO_ASSERT
  assert(v.size()>0);
#endif
  return *std::min_element(v.begin(), v.end());
}

// the + operator

template<class T>
Vector<T> operator+(const Vector<T> &v, T f) 
{
  Vector<T> c(v);
  return c += f;
}

template<class T>
Vector<T> operator+(T f, const Vector<T> &v) 
{
  Vector<T> c(v);
  return c += f;
}

template<class T>
Vector<T> operator+(const Vector<T> &v, const Vector<T> &w) 
{
#ifndef NO_ASSERT
  assert(v.size()==w.size());
#endif
  Vector<T> c(v);
  return c += w;
}

template<class T>
Vector<T> operator-(const Vector<T> &v, T f)
{
  Vector<T> c(v);
  return c -= f;
}

template<class T>
Vector<T> operator-(T f, const Vector<T> &v)
{
  Vector<T> c = -v;
  return c += f;
}

template<class T>
Vector<T> operator-(const Vector<T> &v, const Vector<T> &w)
{
  Vector<T> c(v);
  return c -= w;
}

template<class T>
Vector<T> operator-(const Vector<T> &v)
{
  Vector<T> c(v.size());
  std::transform(v.begin(), v.end(), c.begin(), std::negate<T>());
  return c;
}

template<class T>
T dot(const Vector<T> &v, const Vector<T> &w) 
{
#ifndef NO_ASSERT
  assert(v.size() == w.size());
#endif
  return std::inner_product(v.begin(), v.end(), w.begin(), static_cast<T>(0));
}

template<class T>
Vector<T> operator*(const Vector<T> &v, const Vector<T> &w) 
{
#ifndef NO_ASSERT
  assert(v.size() == w.size());
#endif
  Vector<T> c(v);
  return c*=w;
}

template<class T>
Vector<T> operator*(T f, const Vector<T> &v)
{
  Vector<T> c(v.size(), f);
  return c *= v;
}

template<class T>
Vector<T> operator*(const Vector<T> &v, T f)
{
  Vector<T> c(v.size(), f);
  return c *= v;
}

template<class T>
Vector<T> operator/(const Vector<T> &v, T f) 
{
  Vector<T> c(v);
  return c/=f;
}

template<class T>
Vector<T> operator/(T f, const Vector<T> &v) 
{
  Vector<T> c(v.size());
  std::transform(v.begin(), v.end(), c.begin(), 
  	bind(std::divides<T>(), f, std::placeholders::_1));
  return c;
}

// element by element division

template<class T>
Vector<T> operator/(const Vector<T> &v, const Vector<T> &w) 
{
  Vector<T> c(v);
  return c/=w;
}

// Math functions

template<class T>
Vector<T> sin(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sin ) );
  return b;  
}  

template<class T>
Vector<T> cos(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( cos ) );
  return b;  
}  

template<class T>
Vector<T> tan(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( tan ) );
  return b;  
}  

template<class T>
Vector<T> asin(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( asin ) );
  return b;  
}  

template<class T>
Vector<T> acos(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( acos ) );
  return b;  
}  

template<class T>
Vector<T> atan(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( atan ) );
  return b;  
}  

template<class T>
Vector<T> exp(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( exp ) );
  return b;  
}  

template<class T>
Vector<T> log(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( log ) );
  return b;  
}  

template<class T>
Vector<T> sqrt(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sqrt ) );
  return b;  
}
  
template<class T>
Vector<T> abs(const Vector<T> &v) 
{ 
  Vector<T> b(v.size());  
  std::transform(v.begin(), v.end(), b.begin(), std::abs<T>);
  return b;  
}

template<>
Vector< std::complex<double> > abs(const Vector< std::complex<double> >& v) 
{
  Vector< std::complex<double> > b(v.size());
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::abs(v[i]);
  return b;
}

Vec real(const Vec_cmplx &v) {
  Vec b(v.size());
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::real(v[i]);
  return b;
}

Vec imag(const Vec_cmplx &v) {
  Vec b(v.size());
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::imag(v[i]);
  return b;
}

Vec_cmplx conjug(const Vec_cmplx &v) {
  Vec_cmplx b(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    b[i] = std::complex<double>(std::real(v[i]), -std::imag(v[i]));
  }
  return b;
}

Vec realabs(const Vec_cmplx &v) 
{
  Vec b(v.size());
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::abs(v[i]);
  return b;
}

template<class T>
Vector<T> pow(const Vector<T> &v, T exp) 
{ 
  Vector<T> b(v.size()); 
  std::transform(v.begin(), v.end(), b.begin(), 
  	             bind(Power<T>(), std::placeholders::_1, exp));  
  return b;  
}  

template<class T>
Vector<T> pow(const Vector<T> &v, int exp) 
{ 
  // we've got to find a more efficient way to do this...
  return pow(v, static_cast<T>(exp));
}  

// grid generator
    
template<class T>
Vector<T> space(T xmin, T xmax, size_t n)
{
  Vector<T> a(n);
  for (size_t i = n; i--;)
    a[i] = xmin + (xmax - xmin)*T(i)/T(n - 1);
  return a;
}
    
template<class T>
Vector<T> range(T imin, size_t n) 
{
  Vector<T> a(n);
  for (size_t i = n; i--;)
    a[i] = imin + static_cast<T>(i);
  return a;
}

// max of 2 vectors

template<class T>
Vector<T> max(const Vector<T> &v, const Vector<T> &w) 
{
#ifndef NO_ASSERT
  assert(v.size() == w.size());
#endif
  size_t n = v.size();
  Vector<T> res(n);
  for (size_t i = 0; i < n; ++i)  
    res[i] = v(i) > w(i) ? v(i) : w(i);
  return res;
}

template<class T> 
Vector<T> max(const Vector<T> &v, T f) 
{
  Vector<T> res(v.size());
  std::transform(v.begin(), v.end(), res.begin(), 
  	             bind(Max<T>(), std::placeholders::_1, f));
  return res;
}

template<class T>
Vector<T> max(T f, const Vector<T> &v) 
{
  return max(v, f);
}

// min of 2 arguments

template<class T> 
Vector<T>  min(const Vector<T> &v, const Vector<T> &w) 
{
#ifndef NO_ASSERT
  assert(v.size() == w.size());
#endif
  size_t n = v.size();
  Vector<T> res(n);
  for (size_t i = 0; i < n; ++i)  
    res[i] = v(i) < w(i) ? v(i) : w(i);
  return res;
}

template<class T> Vector<T>  min(const Vector<T> &v, T f) 
{
  Vector<T> res(v.size());
  std::transform(v.begin(), v.end(), res.begin(), 
  	             bind(Min<T>(), std::placeholders::_1, f));
  return res;
}

template<class T>
Vector<T> min(T f, const Vector<T> &v) 
{
  return min(v, f);
}

/*
sum of elements
*/

template<class T>
T sum(const Vector<T> &v) 
{
  return std::accumulate(v.begin(), v.end(), static_cast<T>(0));
}

// concatenate two Vectors

template<class T> Vector<T> cat( const Vector<T> &v1, const Vector<T> &v2) 
{
  size_t n1 = v1.size(); 
  size_t n2 = v2.size();
  Vector<T> res(n1+n2);
  for (size_t i=0; i<n1; ++i) res(i   ) = v1(i);
  for (size_t j=0; j<n2; ++j) res(j+n1) = v2(j);
  return res;
}

// class methods to fill up elements of an existing Vector

template<class T>
void space(Vector<T> &v, T xmin, T xmax) 
{
#ifndef NO_ASSERT
  assert(v.size() >= 2);
#endif
  size_t n = v.size();
  for (size_t i = n; i--;)
    v[i] = xmin + (xmax - xmin)*T(i)/T(n - 1);
}
    
template<class T>
void range(Vector<T> &v, T imin) 
{
  for (size_t i = v.size(); i--;)
    v[i] = T(i) + imin;
}

// return index of min(v)

template<class T>
size_t index_min(Vector<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = min_element(beg, end);
  return static_cast<size_t>(i - beg);
}

// return index of max(v)

template<class T>
size_t index_max(Vector<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = max_element(beg, end);
  return static_cast<size_t>(i - beg);
}

void setReal(Vec_cmplx &v, const Vec& rV) {
  for (size_t i = 0; i < rV.size(); ++i) {
    v[i] = std::complex<double>(rV[i], v[i].imag());
  }
}

void setImag(Vec_cmplx &v, const Vec& iV) {
  for (size_t i = 0; i < iV.size(); ++i) {
    v[i] = std::complex<double>(v[i].real(), iV[i]);
  }
}

void setRealImag(Vec_cmplx &v, const Vec& rV, const Vec& iV) {
  size_t n = rV.size() >= iV.size()? rV.size(): iV.size();
  for (size_t i = 0; i < n; ++i) {
    v[i] = std::complex<double>(rV[i], iV[i]);
  }
}

Vec_cmplx cmplx(const Vec& rV, const Vec& iV) {
  Vec_cmplx res(rV.size());
  setRealImag(res, rV, iV);
  return res;
}

//*******************************************

template <class T>
std::ostream& operator<<(std::ostream& s, const Vector<T>& v)
{
  for_each(v.begin(), v.end(), Print<T>(s));
  return s;
}

//
// template instantiations
//

// double
template class Vector<double>;
template Vector<double> space(double, double, size_t);
template double max(const Vector<double>&);
template Vector<double> max(const Vector<double>&, double);
template Vector<double> max(double, const Vector<double>&);
template Vector<double> max(const Vector<double>&, const Vector<double>&);
template double min(const Vector<double>&);
template Vector<double> min(const Vector<double>&, double);
template Vector<double> min(double, const Vector<double>&);
template Vector<double> min(const Vector<double>&, const Vector<double>&);
template double dot(const Vector<double>&, const Vector<double>&);
template double sum(const Vector<double>&);
template Vector<double> abs(const Vector<double>&);
template Vector<double> pow(const Vector<double>&, double);
template Vector<double> pow(const Vector<double>&, int);
template Vector<double> sqrt(const Vector<double>&);
template Vector<double> exp(const Vector<double>&);
template Vector<double> log(const Vector<double>&);
template Vector<double> sin(const Vector<double>&);
template Vector<double> cos(const Vector<double>&);
template Vector<double> tan(const Vector<double>&);
template Vector<double> asin(const Vector<double>&);
template Vector<double> acos(const Vector<double>&);
template Vector<double> atan(const Vector<double>&);
template std::ostream& operator<<(std::ostream& s, const Vector<double>&);
template Vector<double> cat(const Vector<double>&, const Vector<double>&);
template Vector<double> operator+(const Vector<double>&, double);
template Vector<double> operator+(double, const Vector<double>&);
template Vector<double> operator+(const Vector<double>&, const Vector<double>&);
template Vector<double> operator-(const Vector<double>&, double);
template Vector<double> operator-(double, const Vector<double>&);
template Vector<double> operator-(const Vector<double>&, const Vector<double>&);
template Vector<double> operator-(const Vector<double>&);
template Vector<double> operator*(const Vector<double>&, double);
template Vector<double> operator*(double, const Vector<double>&);
template Vector<double> operator*(const Vector<double>&, const Vector<double>&);
template Vector<double> operator/(const Vector<double>&, const Vector<double>&);
template Vector<double> operator/(const Vector<double>&, double);
template Vector<double> operator/(double, const Vector<double>&);

// size_t
template class Vector<size_t>;
template Vector<size_t> range(size_t, size_t);
template Vector<size_t> space(size_t, size_t, size_t);
template size_t max(const Vector<size_t>&);
template Vector<size_t> max(const Vector<size_t>&, size_t);
template Vector<size_t> max(size_t, const Vector<size_t>&);
template Vector<size_t> max(const Vector<size_t>&, const Vector<size_t>&);
template size_t min(const Vector<size_t>&);
template Vector<size_t> min(const Vector<size_t>&, size_t);
template Vector<size_t> min(size_t, const Vector<size_t>&);
template Vector<size_t> min(const Vector<size_t>&, const Vector<size_t>&);
template size_t dot(const Vector<size_t>&, const Vector<size_t>&);
template size_t sum(const Vector<size_t>&);
template std::ostream& operator<<(std::ostream& s, const Vector<size_t>&);
template Vector<size_t> cat(const Vector<size_t>&, const Vector<size_t>&);
template Vector<size_t> operator+(const Vector<size_t>&, size_t);
template Vector<size_t> operator+(size_t, const Vector<size_t>&);
template Vector<size_t> operator+(const Vector<size_t>&, const Vector<size_t>&);
template Vector<size_t> operator-(const Vector<size_t>&, size_t);
template Vector<size_t> operator-(size_t, const Vector<size_t>&);
template Vector<size_t> operator-(const Vector<size_t>&, const Vector<size_t>&);
template Vector<size_t> operator-(const Vector<size_t>&);
template Vector<size_t> operator*(const Vector<size_t>&, size_t);
template Vector<size_t> operator*(size_t, const Vector<size_t>&);
template Vector<size_t> operator*(const Vector<size_t>&, const Vector<size_t>&);

// int
template class Vector<int>;
template Vector<int> space(int, int, size_t);
template int max(const Vector<int>&);
template Vector<int> max(const Vector<int>&, int);
template Vector<int> max(int, const Vector<int>&);
template Vector<int> max(const Vector<int>&, const Vector<int>&);
template int min(const Vector<int>&);
template Vector<int> min(const Vector<int>&, int);
template Vector<int> min(int, const Vector<int>&);
template Vector<int> min(const Vector<int>&, const Vector<int>&);
template int dot(const Vector<int>&, const Vector<int>&);
template int sum(const Vector<int>&);
template std::ostream& operator<<(std::ostream& s, const Vector<int>&);
template Vector<int> cat(const Vector<int>&, const Vector<int>&);
template Vector<int> operator+(const Vector<int>&, int);
template Vector<int> operator+(int, const Vector<int>&);
template Vector<int> operator+(const Vector<int>&, const Vector<int>&);
template Vector<int> operator-(const Vector<int>&, int);
template Vector<int> operator-(int, const Vector<int>&);
template Vector<int> operator-(const Vector<int>&, const Vector<int>&);
template Vector<int> operator-(const Vector<int>&);
template Vector<int> operator*(const Vector<int>&, int);
template Vector<int> operator*(int, const Vector<int>&);
template Vector<int> operator*(const Vector<int>&, const Vector<int>&);


template class Vector< std::complex<double> >;
template Vector< std::complex<double> > space( std::complex<double> ,  std::complex<double> , size_t);
template std::complex<double>  dot(const Vector< std::complex<double> >&, const Vector< std::complex<double> >&);
template std::complex<double>  sum(const Vector< std::complex<double> >&);
template Vector< std::complex<double> > pow(const Vector< std::complex<double> >&,  std::complex<double> );
template std::ostream& operator<<(std::ostream& s, const Vector< std::complex<double> >&);
template Vector< std::complex<double> > cat(const Vector< std::complex<double> >&, const Vector< std::complex<double> >&);
template Vector< std::complex<double> > operator+(const Vector< std::complex<double> >&,  std::complex<double> );
template Vector< std::complex<double> > operator-(const Vector< std::complex<double> >&,  std::complex<double> );
template Vector< std::complex<double> > operator*(const Vector< std::complex<double> >&,  std::complex<double> );
template Vector< std::complex<double> > operator*( std::complex<double> , const Vector< std::complex<double> >&);
template Vector< std::complex<double> > operator+(const Vector< std::complex<double> >&, const Vector< std::complex<double> >&);
template Vector< std::complex<double> > operator-(const Vector< std::complex<double> >&, const Vector< std::complex<double> >&);
template Vector< std::complex<double> > operator*(const Vector< std::complex<double> >&, const Vector< std::complex<double> >&);
template Vector< std::complex<double> > operator/(const Vector< std::complex<double> >&, const Vector< std::complex<double> >&);
template Vector< std::complex<double> > operator-(const Vector< std::complex<double> >&);


