#include <complex>
#include "mntVec9.h"

/****************************************************************************/
// Definition of the 3-dimensional vector Class
/****************************************************************************/

template<class T>
Vector9<T>::Vector9()
{
}

template<class T>
Vector9<T>::Vector9(T e)
{
  for (size_t i = 0; i < this->size(); ++i)
  {
    (*this)[i] = e;
  }
}

template<class T>
Vector9<T>::Vector9(const Vector9<T> &w)
{
  (*this) = w;
}

template<class T> 
void Vector9<T>::random(void) 
{
  std::transform(this->begin(), this->end(), this->begin(), Random<T>());
}


template<class T>
Vector9<T> &Vector9<T>::operator=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), Setval<T>(f));
  return (*this);
}

template<class T>
Vector9<T> &Vector9<T>::operator+=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        std::bind2nd(std::plus<T>(), f));
  return (*this);
}

template<class T>
Vector9<T> &Vector9<T>::operator+=(const Vector9<T> &w) 
{
  std::transform(this->begin(), this->end(), 
             w.begin(), 
             this->begin(), 
             std::plus<T>() );
  return (*this);
}

template<class T>
Vector9<T> &Vector9<T>::operator-=(const Vector9<T> &w)
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::minus<T>() );
  return (*this);
}

template<class T>
Vector9<T> &Vector9<T>::operator-=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::minus<T>(), f));
  return (*this);
}

template<class T>
Vector9<T> &Vector9<T>::operator*=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::multiplies<T>(), f));
  return *this;  
}

template<class T>
Vector9<T> &Vector9<T>::operator*=(const Vector9<T> &w)
{  
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::multiplies<T>() );
  return (*this);  
}

template<class T>
Vector9<T> &Vector9<T>::operator/=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::divides<T>(), f));
  return (*this);
}

template<class T>
Vector9<T> &Vector9<T>::operator/=(const Vector9<T> &w) 
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::divides<T>() );
  return (*this);
}


// Definition of extra functions supporting the Vector9 class

// find max element

template<class T>
T max(const Vector9<T> &v)
{
  return *std::max_element(v.begin(), v.end());
}

// find min element 

template<class T>
T min(const Vector9<T> &v) 
{
  return *std::min_element(v.begin(), v.end());
}

// the + operator

template<class T>
Vector9<T> operator+(const Vector9<T> &v, T f) 
{
  Vector9<T> c(v);
  return c += f;
}

template<class T>
Vector9<T> operator+(T f, const Vector9<T> &v) 
{
  Vector9<T> c(v);
  return c += f;
}

template<class T>
Vector9<T> operator+(const Vector9<T> &v, const Vector9<T> &w) 
{
  Vector9<T> c(v);
  return c += w;
}

template<class T>
Vector9<T> operator-(const Vector9<T> &v, T f)
{
  Vector9<T> c(v);
  return c -= f;
}

template<class T>
Vector9<T> operator-(T f, const Vector9<T> &v)
{
  Vector9<T> c = -v;
  return c += f;
}

template<class T>
Vector9<T> operator-(const Vector9<T> &v, const Vector9<T> &w)
{
  Vector9<T> c(v);
  return c -= w;
}

template<class T>
Vector9<T> operator-(const Vector9<T> &v)
{
  Vector9<T> c;
  std::transform(v.begin(), v.end(), c.begin(), std::negate<T>());
  return c;
}

template<class T>
T dot(const Vector9<T> &v, const Vector9<T> &w) 
{
  return std::inner_product(v.begin(), v.end(), w.begin(), static_cast<T>(0));
}

template<class T>
Vector9<T> operator*(const Vector9<T> &v, const Vector9<T> &w) 
{
  Vector9<T> c(v);
  return c*=w;
}

template<class T>
Vector9<T> operator*(T f, const Vector9<T> &v)
{
  Vector9<T> c(f);
  return c *= v;
}

template<class T>
Vector9<T> operator*(const Vector9<T> &v, T f)
{
  Vector9<T> c(f);
  return c *= v;
}

template<class T>
Vector9<T> operator/(const Vector9<T> &v, T f) 
{
  Vector9<T> c(v);
  return c/=f;
}

template<class T>
Vector9<T> operator/(T f, const Vector9<T> &v) 
{
  Vector9<T> c;
  std::transform(v.begin(), v.end(), c.begin(), bind1st(std::divides<T>(), f));
  return c;
}

// element by element division

template<class T>
Vector9<T> operator/(const Vector9<T> &v, const Vector9<T> &w) 
{
  Vector9<T> c(v);
  return c/=w;
}

// Math functions

template<class T>
Vector9<T> sin(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sin ) );
  return b;  
}  

template<class T>
Vector9<T> cos(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( cos ) );
  return b;  
}  

template<class T>
Vector9<T> tan(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( tan ) );
  return b;  
}  

template<class T>
Vector9<T> asin(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( asin ) );
  return b;  
}  

template<class T>
Vector9<T> acos(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( acos ) );
  return b;  
}  

template<class T>
Vector9<T> atan(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( atan ) );
  return b;  
}  

template<class T>
Vector9<T> exp(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( exp ) );
  return b;  
}  

template<class T>
Vector9<T> log(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( log ) );
  return b;  
}  

template<class T>
Vector9<T> sqrt(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sqrt ) );
  return b;  
}
  
template<class T>
Vector9<T> abs(const Vector9<T> &v) 
{ 
  Vector9<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), std::abs<T>);
  return b;  
}

template<>
Vector9< std::complex<double> > abs(const Vector9< std::complex<double> >& v) 
{
  Vector9< std::complex<double> > b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::abs(v[i]);
  return b;
}

Vec9 real(const Vec9_cmplx &v) {
  Vec9 b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::real(v[i]);
  return b;
}

Vec9 imag(const Vec9_cmplx &v) {
  Vec9 b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::imag(v[i]);
  return b;
}

Vec9_cmplx conjug(const Vec9_cmplx &v) {
  Vec9_cmplx b;
  for (size_t i = 0; i < v.size(); ++i) {
    b[i] = std::complex<double>(std::real(v[i]), -std::imag(v[i]));
  }
  return b;
}

template<class T>
Vector9<T> pow(const Vector9<T> &v, T exp) 
{ 
  Vector9<T> b; 
  std::transform(v.begin(), v.end(), b.begin(), bind2nd(Power<T>(), exp));  
  return b;  
}  

template<class T>
Vector9<T> pow(const Vector9<T> &v, int exp) 
{ 
  // we've got to find a more efficient way to do this...
  return pow(v, static_cast<T>(exp));
}  

// max of 2 vectors

template<class T>
Vector9<T> max(const Vector9<T> &v, const Vector9<T> &w) 
{
  Vector9<T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) > w(i) ? v(i) : w(i);
  return res;
}

template<class T> 
Vector9<T> max(const Vector9<T> &v, T f) 
{
  Vector9<T> res;
  std::transform(v.begin(), v.end(), res.begin(), bind2nd(Max<T>(), f));
  return res;
}

template<class T>
Vector9<T> max(T f, const Vector9<T> &v) 
{
  return max(v, f);
}

// min of 2 arguments

template<class T> 
Vector9<T>  min(const Vector9<T> &v, const Vector9<T> &w) 
{
  Vector9<T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) < w(i) ? v(i) : w(i);
  return res;
}

template<class T> Vector9<T>  min(const Vector9<T> &v, T f) 
{
  Vector9<T> res;
  std::transform(v.begin(), v.end(), res.begin(), bind2nd(Min<T>(), f));
  return res;
}

template<class T>
Vector9<T> min(T f, const Vector9<T> &v) 
{
  return min(v, f);
}

/*
sum of elements
*/

template<class T>
T sum(const Vector9<T> &v) 
{
  return std::accumulate(v.begin(), v.end(), static_cast<T>(0));
}


// return index of min(v)

template<class T>
size_t index_min(Vector9<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = min_element(beg, end);
  return static_cast<size_t>(i - beg);
}

// return index of max(v)

template<class T>
size_t index_max(Vector9<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = max_element(beg, end);
  return static_cast<size_t>(i - beg);
}

void setReal(Vec9_cmplx &v, const Vec9& rV) {
  for (size_t i = 0; i < rV.size(); ++i) {
    v[i] = std::complex<double>(rV[i], v[i].imag());
  }
}

void setImag(Vec9_cmplx &v, const Vec9& iV) {
  for (size_t i = 0; i < iV.size(); ++i) {
    v[i] = std::complex<double>(v[i].real(), iV[i]);
  }
}

void setRealImag(Vec9_cmplx &v, const Vec9& rV, const Vec9& iV) {
  size_t n = rV.size() >= iV.size()? rV.size(): iV.size();
  for (size_t i = 0; i < n; ++i) {
    v[i] = std::complex<double>(rV[i], iV[i]);
  }
}

Vec9_cmplx cmplx(const Vec9& rV, const Vec9& iV) {
  Vec9_cmplx res;
  setRealImag(res, rV, iV);
  return res;
}

//*******************************************

template <class T>
std::ostream& operator<<(std::ostream& s, const Vector9<T>& v)
{
  for_each(v.begin(), v.end(), Print<T>(s));
  return s;
}

//
// template instantiations
//

// double
template class Vector9<double>;
template double max(const Vector9<double>&);
template Vector9<double> max(const Vector9<double>&, double);
template Vector9<double> max(double, const Vector9<double>&);
template Vector9<double> max(const Vector9<double>&, const Vector9<double>&);
template double min(const Vector9<double>&);
template Vector9<double> min(const Vector9<double>&, double);
template Vector9<double> min(double, const Vector9<double>&);
template Vector9<double> min(const Vector9<double>&, const Vector9<double>&);
template double dot(const Vector9<double>&, const Vector9<double>&);
template double sum(const Vector9<double>&);
template Vector9<double> abs(const Vector9<double>&);
template Vector9<double> pow(const Vector9<double>&, double);
template Vector9<double> pow(const Vector9<double>&, int);
template Vector9<double> sqrt(const Vector9<double>&);
template Vector9<double> exp(const Vector9<double>&);
template Vector9<double> log(const Vector9<double>&);
template Vector9<double> sin(const Vector9<double>&);
template Vector9<double> cos(const Vector9<double>&);
template Vector9<double> tan(const Vector9<double>&);
template Vector9<double> asin(const Vector9<double>&);
template Vector9<double> acos(const Vector9<double>&);
template Vector9<double> atan(const Vector9<double>&);
template std::ostream& operator<<(std::ostream& s, const Vector9<double>&);
template Vector9<double> operator+(const Vector9<double>&, double);
template Vector9<double> operator+(double, const Vector9<double>&);
template Vector9<double> operator+(const Vector9<double>&, const Vector9<double>&);
template Vector9<double> operator-(const Vector9<double>&, double);
template Vector9<double> operator-(double, const Vector9<double>&);
template Vector9<double> operator-(const Vector9<double>&, const Vector9<double>&);
template Vector9<double> operator-(const Vector9<double>&);
template Vector9<double> operator*(const Vector9<double>&, double);
template Vector9<double> operator*(double, const Vector9<double>&);
template Vector9<double> operator*(const Vector9<double>&, const Vector9<double>&);
template Vector9<double> operator/(const Vector9<double>&, const Vector9<double>&);
template Vector9<double> operator/(const Vector9<double>&, double);
template Vector9<double> operator/(double, const Vector9<double>&);

// size_t
template class Vector9<size_t>;
template size_t max(const Vector9<size_t>&);
template Vector9<size_t> max(const Vector9<size_t>&, size_t);
template Vector9<size_t> max(size_t, const Vector9<size_t>&);
template Vector9<size_t> max(const Vector9<size_t>&, const Vector9<size_t>&);
template size_t min(const Vector9<size_t>&);
template Vector9<size_t> min(const Vector9<size_t>&, size_t);
template Vector9<size_t> min(size_t, const Vector9<size_t>&);
template Vector9<size_t> min(const Vector9<size_t>&, const Vector9<size_t>&);
template size_t dot(const Vector9<size_t>&, const Vector9<size_t>&);
template size_t sum(const Vector9<size_t>&);
template std::ostream& operator<<(std::ostream& s, const Vector9<size_t>&);
template Vector9<size_t> operator+(const Vector9<size_t>&, size_t);
template Vector9<size_t> operator+(size_t, const Vector9<size_t>&);
template Vector9<size_t> operator+(const Vector9<size_t>&, const Vector9<size_t>&);
template Vector9<size_t> operator-(const Vector9<size_t>&, size_t);
template Vector9<size_t> operator-(size_t, const Vector9<size_t>&);
template Vector9<size_t> operator-(const Vector9<size_t>&, const Vector9<size_t>&);
template Vector9<size_t> operator-(const Vector9<size_t>&);
template Vector9<size_t> operator*(const Vector9<size_t>&, size_t);
template Vector9<size_t> operator*(size_t, const Vector9<size_t>&);
template Vector9<size_t> operator*(const Vector9<size_t>&, const Vector9<size_t>&);

// int
template class Vector9<int>;
template int max(const Vector9<int>&);
template Vector9<int> max(const Vector9<int>&, int);
template Vector9<int> max(int, const Vector9<int>&);
template Vector9<int> max(const Vector9<int>&, const Vector9<int>&);
template int min(const Vector9<int>&);
template Vector9<int> min(const Vector9<int>&, int);
template Vector9<int> min(int, const Vector9<int>&);
template Vector9<int> min(const Vector9<int>&, const Vector9<int>&);
template int dot(const Vector9<int>&, const Vector9<int>&);
template int sum(const Vector9<int>&);
template std::ostream& operator<<(std::ostream& s, const Vector9<int>&);
template Vector9<int> operator+(const Vector9<int>&, int);
template Vector9<int> operator+(int, const Vector9<int>&);
template Vector9<int> operator+(const Vector9<int>&, const Vector9<int>&);
template Vector9<int> operator-(const Vector9<int>&, int);
template Vector9<int> operator-(int, const Vector9<int>&);
template Vector9<int> operator-(const Vector9<int>&, const Vector9<int>&);
template Vector9<int> operator-(const Vector9<int>&);
template Vector9<int> operator*(const Vector9<int>&, int);
template Vector9<int> operator*(int, const Vector9<int>&);
template Vector9<int> operator*(const Vector9<int>&, const Vector9<int>&);


template class Vector9< std::complex<double> >;
template std::complex<double>  dot(const Vector9< std::complex<double> >&, const Vector9< std::complex<double> >&);
template std::complex<double>  sum(const Vector9< std::complex<double> >&);
template Vector9< std::complex<double> > pow(const Vector9< std::complex<double> >&,  std::complex<double> );
template std::ostream& operator<<(std::ostream& s, const Vector9< std::complex<double> >&);
template Vector9< std::complex<double> > operator+(const Vector9< std::complex<double> >&,  std::complex<double> );
template Vector9< std::complex<double> > operator-(const Vector9< std::complex<double> >&,  std::complex<double> );
template Vector9< std::complex<double> > operator*(const Vector9< std::complex<double> >&,  std::complex<double> );
template Vector9< std::complex<double> > operator*( std::complex<double> , const Vector9< std::complex<double> >&);
template Vector9< std::complex<double> > operator+(const Vector9< std::complex<double> >&, const Vector9< std::complex<double> >&);
template Vector9< std::complex<double> > operator-(const Vector9< std::complex<double> >&, const Vector9< std::complex<double> >&);
template Vector9< std::complex<double> > operator*(const Vector9< std::complex<double> >&, const Vector9< std::complex<double> >&);
template Vector9< std::complex<double> > operator/(const Vector9< std::complex<double> >&, const Vector9< std::complex<double> >&);
template Vector9< std::complex<double> > operator-(const Vector9< std::complex<double> >&);


