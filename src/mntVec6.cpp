#include <complex>
#include "mntVec6.h"

/****************************************************************************/
// Definition of the 3-dimensional vector Class
/****************************************************************************/

template<class T>
Vector6<T>::Vector6()
{
}

template<class T>
Vector6<T>::Vector6(T e)
{
  for (size_t i = 0; i < this->size(); ++i)
  {
    (*this)[i] = e;
  }
}

template<class T>
Vector6<T>::Vector6(T* ptr)
{
  for (size_t i = 0; i < this->size(); ++i)
  {
    (*this)[i] = *(ptr + i);
  }
}

template<class T>
Vector6<T>::Vector6(const Vector6<T> &w)
{
  (*this) = w;
}

template<class T> 
void Vector6<T>::random(void) 
{
  std::transform(this->begin(), this->end(), this->begin(), Random<T>());
}


template<class T>
Vector6<T> &Vector6<T>::operator=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), Setval<T>(f));
  return (*this);
}

template<class T>
Vector6<T> &Vector6<T>::operator+=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        std::bind2nd(std::plus<T>(), f));
  return (*this);
}

template<class T>
Vector6<T> &Vector6<T>::operator+=(const Vector6<T> &w) 
{
  std::transform(this->begin(), this->end(), 
             w.begin(), 
             this->begin(), 
             std::plus<T>() );
  return (*this);
}

template<class T>
Vector6<T> &Vector6<T>::operator-=(const Vector6<T> &w)
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::minus<T>() );
  return (*this);
}

template<class T>
Vector6<T> &Vector6<T>::operator-=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::minus<T>(), f));
  return (*this);
}

template<class T>
Vector6<T> &Vector6<T>::operator*=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::multiplies<T>(), f));
  return *this;  
}

template<class T>
Vector6<T> &Vector6<T>::operator*=(const Vector6<T> &w)
{  
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::multiplies<T>() );
  return (*this);  
}

template<class T>
Vector6<T> &Vector6<T>::operator/=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::divides<T>(), f));
  return (*this);
}

template<class T>
Vector6<T> &Vector6<T>::operator/=(const Vector6<T> &w) 
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::divides<T>() );
  return (*this);
}


// Definition of extra functions supporting the Vector6 class

// find max element

template<class T>
T max(const Vector6<T> &v)
{
  return *std::max_element(v.begin(), v.end());
}

// find min element 

template<class T>
T min(const Vector6<T> &v) 
{
  return *std::min_element(v.begin(), v.end());
}

// the + operator

template<class T>
Vector6<T> operator+(const Vector6<T> &v, T f) 
{
  Vector6<T> c(v);
  return c += f;
}

template<class T>
Vector6<T> operator+(T f, const Vector6<T> &v) 
{
  Vector6<T> c(v);
  return c += f;
}

template<class T>
Vector6<T> operator+(const Vector6<T> &v, const Vector6<T> &w) 
{
  Vector6<T> c(v);
  return c += w;
}

template<class T>
Vector6<T> operator-(const Vector6<T> &v, T f)
{
  Vector6<T> c(v);
  return c -= f;
}

template<class T>
Vector6<T> operator-(T f, const Vector6<T> &v)
{
  Vector6<T> c = -v;
  return c += f;
}

template<class T>
Vector6<T> operator-(const Vector6<T> &v, const Vector6<T> &w)
{
  Vector6<T> c(v);
  return c -= w;
}

template<class T>
Vector6<T> operator-(const Vector6<T> &v)
{
  Vector6<T> c;
  std::transform(v.begin(), v.end(), c.begin(), std::negate<T>());
  return c;
}

template<class T>
T dot(const Vector6<T> &v, const Vector6<T> &w) 
{
  return std::inner_product(v.begin(), v.end(), w.begin(), static_cast<T>(0));
}

template<class T>
Vector6<T> operator*(const Vector6<T> &v, const Vector6<T> &w) 
{
  Vector6<T> c(v);
  return c*=w;
}

template<class T>
Vector6<T> operator*(T f, const Vector6<T> &v)
{
  Vector6<T> c(f);
  return c *= v;
}

template<class T>
Vector6<T> operator*(const Vector6<T> &v, T f)
{
  Vector6<T> c(f);
  return c *= v;
}

template<class T>
Vector6<T> operator/(const Vector6<T> &v, T f) 
{
  Vector6<T> c(v);
  return c/=f;
}

template<class T>
Vector6<T> operator/(T f, const Vector6<T> &v) 
{
  Vector6<T> c;
  std::transform(v.begin(), v.end(), c.begin(), bind1st(std::divides<T>(), f));
  return c;
}

// element by element division

template<class T>
Vector6<T> operator/(const Vector6<T> &v, const Vector6<T> &w) 
{
  Vector6<T> c(v);
  return c/=w;
}

// Math functions

template<class T>
Vector6<T> sin(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sin ) );
  return b;  
}  

template<class T>
Vector6<T> cos(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( cos ) );
  return b;  
}  

template<class T>
Vector6<T> tan(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( tan ) );
  return b;  
}  

template<class T>
Vector6<T> asin(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( asin ) );
  return b;  
}  

template<class T>
Vector6<T> acos(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( acos ) );
  return b;  
}  

template<class T>
Vector6<T> atan(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( atan ) );
  return b;  
}  

template<class T>
Vector6<T> exp(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( exp ) );
  return b;  
}  

template<class T>
Vector6<T> log(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( log ) );
  return b;  
}  

template<class T>
Vector6<T> sqrt(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sqrt ) );
  return b;  
}
  
template<class T>
Vector6<T> abs(const Vector6<T> &v) 
{ 
  Vector6<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), std::abs<T>);
  return b;  
}

template<>
Vector6< std::complex<double> > abs(const Vector6< std::complex<double> >& v) 
{
  Vector6< std::complex<double> > b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::abs(v[i]);
  return b;
}

Vec6 real(const Vec6_cmplx &v) {
  Vec6 b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::real(v[i]);
  return b;
}

Vec6 imag(const Vec6_cmplx &v) {
  Vec6 b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::imag(v[i]);
  return b;
}

Vec6_cmplx conjug(const Vec6_cmplx &v) {
  Vec6_cmplx b;
  for (size_t i = 0; i < v.size(); ++i) {
    b[i] = std::complex<double>(std::real(v[i]), -std::imag(v[i]));
  }
  return b;
}

template<class T>
Vector6<T> pow(const Vector6<T> &v, T exp) 
{ 
  Vector6<T> b; 
  std::transform(v.begin(), v.end(), b.begin(), bind2nd(Power<T>(), exp));  
  return b;  
}  

template<class T>
Vector6<T> pow(const Vector6<T> &v, int exp) 
{ 
  // we've got to find a more efficient way to do this...
  return pow(v, static_cast<T>(exp));
}  

// max of 2 vectors

template<class T>
Vector6<T> max(const Vector6<T> &v, const Vector6<T> &w) 
{
  Vector6<T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) > w(i) ? v(i) : w(i);
  return res;
}

template<class T> 
Vector6<T> max(const Vector6<T> &v, T f) 
{
  Vector6<T> res;
  std::transform(v.begin(), v.end(), res.begin(), bind2nd(Max<T>(), f));
  return res;
}

template<class T>
Vector6<T> max(T f, const Vector6<T> &v) 
{
  return max(v, f);
}

// min of 2 arguments

template<class T> 
Vector6<T>  min(const Vector6<T> &v, const Vector6<T> &w) 
{
  Vector6<T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) < w(i) ? v(i) : w(i);
  return res;
}

template<class T> Vector6<T>  min(const Vector6<T> &v, T f) 
{
  Vector6<T> res;
  std::transform(v.begin(), v.end(), res.begin(), bind2nd(Min<T>(), f));
  return res;
}

template<class T>
Vector6<T> min(T f, const Vector6<T> &v) 
{
  return min(v, f);
}

/*
sum of elements
*/

template<class T>
T sum(const Vector6<T> &v) 
{
  return std::accumulate(v.begin(), v.end(), static_cast<T>(0));
}


// return index of min(v)

template<class T>
size_t index_min(Vector6<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = min_element(beg, end);
  return static_cast<size_t>(i - beg);
}

// return index of max(v)

template<class T>
size_t index_max(Vector6<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = max_element(beg, end);
  return static_cast<size_t>(i - beg);
}

void setReal(Vec6_cmplx &v, const Vec6& rV) {
  for (size_t i = 0; i < rV.size(); ++i) {
    v[i] = std::complex<double>(rV[i], v[i].imag());
  }
}

void setImag(Vec6_cmplx &v, const Vec6& iV) {
  for (size_t i = 0; i < iV.size(); ++i) {
    v[i] = std::complex<double>(v[i].real(), iV[i]);
  }
}

void setRealImag(Vec6_cmplx &v, const Vec6& rV, const Vec6& iV) {
  size_t n = rV.size() >= iV.size()? rV.size(): iV.size();
  for (size_t i = 0; i < n; ++i) {
    v[i] = std::complex<double>(rV[i], iV[i]);
  }
}

Vec6_cmplx cmplx(const Vec6& rV, const Vec6& iV) {
  Vec6_cmplx res;
  setRealImag(res, rV, iV);
  return res;
}

//*******************************************

template <class T>
std::ostream& operator<<(std::ostream& s, const Vector6<T>& v)
{
  for_each(v.begin(), v.end(), Print<T>(s));
  return s;
}

//
// template instantiations
//

// double
template class Vector6<double>;
template double max(const Vector6<double>&);
template Vector6<double> max(const Vector6<double>&, double);
template Vector6<double> max(double, const Vector6<double>&);
template Vector6<double> max(const Vector6<double>&, const Vector6<double>&);
template double min(const Vector6<double>&);
template Vector6<double> min(const Vector6<double>&, double);
template Vector6<double> min(double, const Vector6<double>&);
template Vector6<double> min(const Vector6<double>&, const Vector6<double>&);
template double dot(const Vector6<double>&, const Vector6<double>&);
template double sum(const Vector6<double>&);
template Vector6<double> abs(const Vector6<double>&);
template Vector6<double> pow(const Vector6<double>&, double);
template Vector6<double> pow(const Vector6<double>&, int);
template Vector6<double> sqrt(const Vector6<double>&);
template Vector6<double> exp(const Vector6<double>&);
template Vector6<double> log(const Vector6<double>&);
template Vector6<double> sin(const Vector6<double>&);
template Vector6<double> cos(const Vector6<double>&);
template Vector6<double> tan(const Vector6<double>&);
template Vector6<double> asin(const Vector6<double>&);
template Vector6<double> acos(const Vector6<double>&);
template Vector6<double> atan(const Vector6<double>&);
template std::ostream& operator<<(std::ostream& s, const Vector6<double>&);
template Vector6<double> operator+(const Vector6<double>&, double);
template Vector6<double> operator+(double, const Vector6<double>&);
template Vector6<double> operator+(const Vector6<double>&, const Vector6<double>&);
template Vector6<double> operator-(const Vector6<double>&, double);
template Vector6<double> operator-(double, const Vector6<double>&);
template Vector6<double> operator-(const Vector6<double>&, const Vector6<double>&);
template Vector6<double> operator-(const Vector6<double>&);
template Vector6<double> operator*(const Vector6<double>&, double);
template Vector6<double> operator*(double, const Vector6<double>&);
template Vector6<double> operator*(const Vector6<double>&, const Vector6<double>&);
template Vector6<double> operator/(const Vector6<double>&, const Vector6<double>&);
template Vector6<double> operator/(const Vector6<double>&, double);
template Vector6<double> operator/(double, const Vector6<double>&);

// size_t
template class Vector6<size_t>;
template size_t max(const Vector6<size_t>&);
template Vector6<size_t> max(const Vector6<size_t>&, size_t);
template Vector6<size_t> max(size_t, const Vector6<size_t>&);
template Vector6<size_t> max(const Vector6<size_t>&, const Vector6<size_t>&);
template size_t min(const Vector6<size_t>&);
template Vector6<size_t> min(const Vector6<size_t>&, size_t);
template Vector6<size_t> min(size_t, const Vector6<size_t>&);
template Vector6<size_t> min(const Vector6<size_t>&, const Vector6<size_t>&);
template size_t dot(const Vector6<size_t>&, const Vector6<size_t>&);
template size_t sum(const Vector6<size_t>&);
template std::ostream& operator<<(std::ostream& s, const Vector6<size_t>&);
template Vector6<size_t> operator+(const Vector6<size_t>&, size_t);
template Vector6<size_t> operator+(size_t, const Vector6<size_t>&);
template Vector6<size_t> operator+(const Vector6<size_t>&, const Vector6<size_t>&);
template Vector6<size_t> operator-(const Vector6<size_t>&, size_t);
template Vector6<size_t> operator-(size_t, const Vector6<size_t>&);
template Vector6<size_t> operator-(const Vector6<size_t>&, const Vector6<size_t>&);
template Vector6<size_t> operator-(const Vector6<size_t>&);
template Vector6<size_t> operator*(const Vector6<size_t>&, size_t);
template Vector6<size_t> operator*(size_t, const Vector6<size_t>&);
template Vector6<size_t> operator*(const Vector6<size_t>&, const Vector6<size_t>&);

// int
template class Vector6<int>;
template int max(const Vector6<int>&);
template Vector6<int> max(const Vector6<int>&, int);
template Vector6<int> max(int, const Vector6<int>&);
template Vector6<int> max(const Vector6<int>&, const Vector6<int>&);
template int min(const Vector6<int>&);
template Vector6<int> min(const Vector6<int>&, int);
template Vector6<int> min(int, const Vector6<int>&);
template Vector6<int> min(const Vector6<int>&, const Vector6<int>&);
template int dot(const Vector6<int>&, const Vector6<int>&);
template int sum(const Vector6<int>&);
template std::ostream& operator<<(std::ostream& s, const Vector6<int>&);
template Vector6<int> operator+(const Vector6<int>&, int);
template Vector6<int> operator+(int, const Vector6<int>&);
template Vector6<int> operator+(const Vector6<int>&, const Vector6<int>&);
template Vector6<int> operator-(const Vector6<int>&, int);
template Vector6<int> operator-(int, const Vector6<int>&);
template Vector6<int> operator-(const Vector6<int>&, const Vector6<int>&);
template Vector6<int> operator-(const Vector6<int>&);
template Vector6<int> operator*(const Vector6<int>&, int);
template Vector6<int> operator*(int, const Vector6<int>&);
template Vector6<int> operator*(const Vector6<int>&, const Vector6<int>&);


template class Vector6< std::complex<double> >;
template std::complex<double>  dot(const Vector6< std::complex<double> >&, const Vector6< std::complex<double> >&);
template std::complex<double>  sum(const Vector6< std::complex<double> >&);
template Vector6< std::complex<double> > pow(const Vector6< std::complex<double> >&,  std::complex<double> );
template std::ostream& operator<<(std::ostream& s, const Vector6< std::complex<double> >&);
template Vector6< std::complex<double> > operator+(const Vector6< std::complex<double> >&,  std::complex<double> );
template Vector6< std::complex<double> > operator-(const Vector6< std::complex<double> >&,  std::complex<double> );
template Vector6< std::complex<double> > operator*(const Vector6< std::complex<double> >&,  std::complex<double> );
template Vector6< std::complex<double> > operator*( std::complex<double> , const Vector6< std::complex<double> >&);
template Vector6< std::complex<double> > operator+(const Vector6< std::complex<double> >&, const Vector6< std::complex<double> >&);
template Vector6< std::complex<double> > operator-(const Vector6< std::complex<double> >&, const Vector6< std::complex<double> >&);
template Vector6< std::complex<double> > operator*(const Vector6< std::complex<double> >&, const Vector6< std::complex<double> >&);
template Vector6< std::complex<double> > operator/(const Vector6< std::complex<double> >&, const Vector6< std::complex<double> >&);
template Vector6< std::complex<double> > operator-(const Vector6< std::complex<double> >&);


