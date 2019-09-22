#include <complex>
#include "mntVec2.h"

/****************************************************************************/
// Definition of the 3-dimensional vector Class
/****************************************************************************/

template<class T>
Vector2<T>::Vector2()
{
}

template<class T>
Vector2<T>::Vector2(T e)
{
  for (size_t i = 0; i < this->size(); ++i)
  {
    (*this)[i] = e;
  }
}

template<class T>
Vector2<T>::Vector2(T* ptr)
{
  for (size_t i = 0; i < this->size(); ++i)
  {
    (*this)[i] = *(ptr + i);
  }
}

template<class T>
Vector2<T>::Vector2(const Vector2<T> &w)
{
  (*this) = w;
}

template<class T> 
void Vector2<T>::random(void) 
{
  std::transform(this->begin(), this->end(), this->begin(), Random<T>());
}


template<class T>
Vector2<T> &Vector2<T>::operator=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), Setval<T>(f));
  return (*this);
}

template<class T>
Vector2<T> &Vector2<T>::operator+=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        std::bind2nd(std::plus<T>(), f));
  return (*this);
}

template<class T>
Vector2<T> &Vector2<T>::operator+=(const Vector2<T> &w) 
{
  std::transform(this->begin(), this->end(), 
             w.begin(), 
             this->begin(), 
             std::plus<T>() );
  return (*this);
}

template<class T>
Vector2<T> &Vector2<T>::operator-=(const Vector2<T> &w)
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::minus<T>() );
  return (*this);
}

template<class T>
Vector2<T> &Vector2<T>::operator-=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::minus<T>(), f));
  return (*this);
}

template<class T>
Vector2<T> &Vector2<T>::operator*=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::multiplies<T>(), f));
  return *this;  
}

template<class T>
Vector2<T> &Vector2<T>::operator*=(const Vector2<T> &w)
{  
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::multiplies<T>() );
  return (*this);  
}

template<class T>
Vector2<T> &Vector2<T>::operator/=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        bind2nd(std::divides<T>(), f));
  return (*this);
}

template<class T>
Vector2<T> &Vector2<T>::operator/=(const Vector2<T> &w) 
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::divides<T>() );
  return (*this);
}


// Definition of extra functions supporting the Vector2 class

// find max element

template<class T>
T max(const Vector2<T> &v)
{
  return *std::max_element(v.begin(), v.end());
}

// find min element 

template<class T>
T min(const Vector2<T> &v) 
{
  return *std::min_element(v.begin(), v.end());
}

// the + operator

template<class T>
Vector2<T> operator+(const Vector2<T> &v, T f) 
{
  Vector2<T> c(v);
  return c += f;
}

template<class T>
Vector2<T> operator+(T f, const Vector2<T> &v) 
{
  Vector2<T> c(v);
  return c += f;
}

template<class T>
Vector2<T> operator+(const Vector2<T> &v, const Vector2<T> &w) 
{
  Vector2<T> c(v);
  return c += w;
}

template<class T>
Vector2<T> operator-(const Vector2<T> &v, T f)
{
  Vector2<T> c(v);
  return c -= f;
}

template<class T>
Vector2<T> operator-(T f, const Vector2<T> &v)
{
  Vector2<T> c = -v;
  return c += f;
}

template<class T>
Vector2<T> operator-(const Vector2<T> &v, const Vector2<T> &w)
{
  Vector2<T> c(v);
  return c -= w;
}

template<class T>
Vector2<T> operator-(const Vector2<T> &v)
{
  Vector2<T> c;
  std::transform(v.begin(), v.end(), c.begin(), std::negate<T>());
  return c;
}

template<class T>
T dot(const Vector2<T> &v, const Vector2<T> &w) 
{
  return std::inner_product(v.begin(), v.end(), w.begin(), static_cast<T>(0));
}

template<class T>
Vector2<T> operator*(const Vector2<T> &v, const Vector2<T> &w) 
{
  Vector2<T> c(v);
  return c*=w;
}

template<class T>
Vector2<T> operator*(T f, const Vector2<T> &v)
{
  Vector2<T> c(f);
  return c *= v;
}

template<class T>
Vector2<T> operator*(const Vector2<T> &v, T f)
{
  Vector2<T> c(f);
  return c *= v;
}

template<class T>
Vector2<T> operator/(const Vector2<T> &v, T f) 
{
  Vector2<T> c(v);
  return c/=f;
}

template<class T>
Vector2<T> operator/(T f, const Vector2<T> &v) 
{
  Vector2<T> c;
  std::transform(v.begin(), v.end(), c.begin(), bind1st(std::divides<T>(), f));
  return c;
}

// element by element division

template<class T>
Vector2<T> operator/(const Vector2<T> &v, const Vector2<T> &w) 
{
  Vector2<T> c(v);
  return c/=w;
}

// Math functions

template<class T>
Vector2<T> sin(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sin ) );
  return b;  
}  

template<class T>
Vector2<T> cos(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( cos ) );
  return b;  
}  

template<class T>
Vector2<T> tan(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( tan ) );
  return b;  
}  

template<class T>
Vector2<T> asin(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( asin ) );
  return b;  
}  

template<class T>
Vector2<T> acos(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( acos ) );
  return b;  
}  

template<class T>
Vector2<T> atan(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( atan ) );
  return b;  
}  

template<class T>
Vector2<T> exp(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( exp ) );
  return b;  
}  

template<class T>
Vector2<T> log(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( log ) );
  return b;  
}  

template<class T>
Vector2<T> sqrt(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sqrt ) );
  return b;  
}
  
template<class T>
Vector2<T> abs(const Vector2<T> &v) 
{ 
  Vector2<T> b;  
  std::transform(v.begin(), v.end(), b.begin(), std::abs<T>);
  return b;  
}

template<>
Vector2< std::complex<double> > abs(const Vector2< std::complex<double> >& v) 
{
  Vector2< std::complex<double> > b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::abs(v[i]);
  return b;
}

Vec2 real(const Vec2_cmplx &v) {
  Vec2 b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::real(v[i]);
  return b;
}

Vec2 imag(const Vec2_cmplx &v) {
  Vec2 b;
  for (size_t i = 0; i < v.size(); ++i)
    b[i] = std::imag(v[i]);
  return b;
}

Vec2_cmplx conjug(const Vec2_cmplx &v) {
  Vec2_cmplx b;
  for (size_t i = 0; i < v.size(); ++i) {
    b[i] = std::complex<double>(std::real(v[i]), -std::imag(v[i]));
  }
  return b;
}

template<class T>
Vector2<T> pow(const Vector2<T> &v, T exp) 
{ 
  Vector2<T> b; 
  std::transform(v.begin(), v.end(), b.begin(), bind2nd(Power<T>(), exp));  
  return b;  
}  

template<class T>
Vector2<T> pow(const Vector2<T> &v, int exp) 
{ 
  // we've got to find a more efficient way to do this...
  return pow(v, static_cast<T>(exp));
}  

// max of 2 vectors

template<class T>
Vector2<T> max(const Vector2<T> &v, const Vector2<T> &w) 
{
  Vector2<T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) > w(i) ? v(i) : w(i);
  return res;
}

template<class T> 
Vector2<T> max(const Vector2<T> &v, T f) 
{
  Vector2<T> res;
  std::transform(v.begin(), v.end(), res.begin(), bind2nd(Max<T>(), f));
  return res;
}

template<class T>
Vector2<T> max(T f, const Vector2<T> &v) 
{
  return max(v, f);
}

// min of 2 arguments

template<class T> 
Vector2<T>  min(const Vector2<T> &v, const Vector2<T> &w) 
{
  Vector2<T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) < w(i) ? v(i) : w(i);
  return res;
}

template<class T> Vector2<T>  min(const Vector2<T> &v, T f) 
{
  Vector2<T> res;
  std::transform(v.begin(), v.end(), res.begin(), bind2nd(Min<T>(), f));
  return res;
}

template<class T>
Vector2<T> min(T f, const Vector2<T> &v) 
{
  return min(v, f);
}

/*
sum of elements
*/

template<class T>
T sum(const Vector2<T> &v) 
{
  return std::accumulate(v.begin(), v.end(), static_cast<T>(0));
}


// return index of min(v)

template<class T>
size_t index_min(Vector2<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = min_element(beg, end);
  return static_cast<size_t>(i - beg);
}

// return index of max(v)

template<class T>
size_t index_max(Vector2<T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = max_element(beg, end);
  return static_cast<size_t>(i - beg);
}

void setReal(Vec2_cmplx &v, const Vec2& rV) {
  for (size_t i = 0; i < rV.size(); ++i) {
    v[i] = std::complex<double>(rV[i], v[i].imag());
  }
}

void setImag(Vec2_cmplx &v, const Vec2& iV) {
  for (size_t i = 0; i < iV.size(); ++i) {
    v[i] = std::complex<double>(v[i].real(), iV[i]);
  }
}

void setRealImag(Vec2_cmplx &v, const Vec2& rV, const Vec2& iV) {
  size_t n = rV.size() >= iV.size()? rV.size(): iV.size();
  for (size_t i = 0; i < n; ++i) {
    v[i] = std::complex<double>(rV[i], iV[i]);
  }
}

Vec2_cmplx cmplx(const Vec2& rV, const Vec2& iV) {
  Vec2_cmplx res;
  setRealImag(res, rV, iV);
  return res;
}

//*******************************************

template <class T>
std::ostream& operator<<(std::ostream& s, const Vector2<T>& v)
{
  for_each(v.begin(), v.end(), Print<T>(s));
  return s;
}

//
// template instantiations
//

// double
template class Vector2<double>;
template double max(const Vector2<double>&);
template Vector2<double> max(const Vector2<double>&, double);
template Vector2<double> max(double, const Vector2<double>&);
template Vector2<double> max(const Vector2<double>&, const Vector2<double>&);
template double min(const Vector2<double>&);
template Vector2<double> min(const Vector2<double>&, double);
template Vector2<double> min(double, const Vector2<double>&);
template Vector2<double> min(const Vector2<double>&, const Vector2<double>&);
template double dot(const Vector2<double>&, const Vector2<double>&);
template double sum(const Vector2<double>&);
template Vector2<double> abs(const Vector2<double>&);
template Vector2<double> pow(const Vector2<double>&, double);
template Vector2<double> pow(const Vector2<double>&, int);
template Vector2<double> sqrt(const Vector2<double>&);
template Vector2<double> exp(const Vector2<double>&);
template Vector2<double> log(const Vector2<double>&);
template Vector2<double> sin(const Vector2<double>&);
template Vector2<double> cos(const Vector2<double>&);
template Vector2<double> tan(const Vector2<double>&);
template Vector2<double> asin(const Vector2<double>&);
template Vector2<double> acos(const Vector2<double>&);
template Vector2<double> atan(const Vector2<double>&);
template std::ostream& operator<<(std::ostream& s, const Vector2<double>&);
template Vector2<double> operator+(const Vector2<double>&, double);
template Vector2<double> operator+(double, const Vector2<double>&);
template Vector2<double> operator+(const Vector2<double>&, const Vector2<double>&);
template Vector2<double> operator-(const Vector2<double>&, double);
template Vector2<double> operator-(double, const Vector2<double>&);
template Vector2<double> operator-(const Vector2<double>&, const Vector2<double>&);
template Vector2<double> operator-(const Vector2<double>&);
template Vector2<double> operator*(const Vector2<double>&, double);
template Vector2<double> operator*(double, const Vector2<double>&);
template Vector2<double> operator*(const Vector2<double>&, const Vector2<double>&);
template Vector2<double> operator/(const Vector2<double>&, const Vector2<double>&);
template Vector2<double> operator/(const Vector2<double>&, double);
template Vector2<double> operator/(double, const Vector2<double>&);

// size_t
template class Vector2<size_t>;
template size_t max(const Vector2<size_t>&);
template Vector2<size_t> max(const Vector2<size_t>&, size_t);
template Vector2<size_t> max(size_t, const Vector2<size_t>&);
template Vector2<size_t> max(const Vector2<size_t>&, const Vector2<size_t>&);
template size_t min(const Vector2<size_t>&);
template Vector2<size_t> min(const Vector2<size_t>&, size_t);
template Vector2<size_t> min(size_t, const Vector2<size_t>&);
template Vector2<size_t> min(const Vector2<size_t>&, const Vector2<size_t>&);
template size_t dot(const Vector2<size_t>&, const Vector2<size_t>&);
template size_t sum(const Vector2<size_t>&);
template std::ostream& operator<<(std::ostream& s, const Vector2<size_t>&);
template Vector2<size_t> operator+(const Vector2<size_t>&, size_t);
template Vector2<size_t> operator+(size_t, const Vector2<size_t>&);
template Vector2<size_t> operator+(const Vector2<size_t>&, const Vector2<size_t>&);
template Vector2<size_t> operator-(const Vector2<size_t>&, size_t);
template Vector2<size_t> operator-(size_t, const Vector2<size_t>&);
template Vector2<size_t> operator-(const Vector2<size_t>&, const Vector2<size_t>&);
template Vector2<size_t> operator-(const Vector2<size_t>&);
template Vector2<size_t> operator*(const Vector2<size_t>&, size_t);
template Vector2<size_t> operator*(size_t, const Vector2<size_t>&);
template Vector2<size_t> operator*(const Vector2<size_t>&, const Vector2<size_t>&);

// int
template class Vector2<int>;
template int max(const Vector2<int>&);
template Vector2<int> max(const Vector2<int>&, int);
template Vector2<int> max(int, const Vector2<int>&);
template Vector2<int> max(const Vector2<int>&, const Vector2<int>&);
template int min(const Vector2<int>&);
template Vector2<int> min(const Vector2<int>&, int);
template Vector2<int> min(int, const Vector2<int>&);
template Vector2<int> min(const Vector2<int>&, const Vector2<int>&);
template int dot(const Vector2<int>&, const Vector2<int>&);
template int sum(const Vector2<int>&);
template std::ostream& operator<<(std::ostream& s, const Vector2<int>&);
template Vector2<int> operator+(const Vector2<int>&, int);
template Vector2<int> operator+(int, const Vector2<int>&);
template Vector2<int> operator+(const Vector2<int>&, const Vector2<int>&);
template Vector2<int> operator-(const Vector2<int>&, int);
template Vector2<int> operator-(int, const Vector2<int>&);
template Vector2<int> operator-(const Vector2<int>&, const Vector2<int>&);
template Vector2<int> operator-(const Vector2<int>&);
template Vector2<int> operator*(const Vector2<int>&, int);
template Vector2<int> operator*(int, const Vector2<int>&);
template Vector2<int> operator*(const Vector2<int>&, const Vector2<int>&);


template class Vector2< std::complex<double> >;
template std::complex<double>  dot(const Vector2< std::complex<double> >&, const Vector2< std::complex<double> >&);
template std::complex<double>  sum(const Vector2< std::complex<double> >&);
template Vector2< std::complex<double> > pow(const Vector2< std::complex<double> >&,  std::complex<double> );
template std::ostream& operator<<(std::ostream& s, const Vector2< std::complex<double> >&);
template Vector2< std::complex<double> > operator+(const Vector2< std::complex<double> >&,  std::complex<double> );
template Vector2< std::complex<double> > operator-(const Vector2< std::complex<double> >&,  std::complex<double> );
template Vector2< std::complex<double> > operator*(const Vector2< std::complex<double> >&,  std::complex<double> );
template Vector2< std::complex<double> > operator*( std::complex<double> , const Vector2< std::complex<double> >&);
template Vector2< std::complex<double> > operator+(const Vector2< std::complex<double> >&, const Vector2< std::complex<double> >&);
template Vector2< std::complex<double> > operator-(const Vector2< std::complex<double> >&, const Vector2< std::complex<double> >&);
template Vector2< std::complex<double> > operator*(const Vector2< std::complex<double> >&, const Vector2< std::complex<double> >&);
template Vector2< std::complex<double> > operator/(const Vector2< std::complex<double> >&, const Vector2< std::complex<double> >&);
template Vector2< std::complex<double> > operator-(const Vector2< std::complex<double> >&);


