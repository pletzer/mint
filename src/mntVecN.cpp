#include "mntVecN.h"
#include <cmath>

template<size_t N, class T>
VecN<N, T>::VecN()
{
}

template<size_t N, class T>
VecN<N, T>::VecN(T e)
{
  for (size_t i = 0; i < this->size(); ++i)
  {
    (*this)[i] = e;
  }
}

template<size_t N, class T>
VecN<N, T>::VecN(const T* ptr)
{
  for (size_t i = 0; i < this->size(); ++i)
  {
    (*this)[i] = *(ptr + i);
  }
}

template<size_t N, class T>
VecN<N, T>::VecN(const VecN<N, T> &w)
{
  (*this) = w;
}

template<size_t N, class T> 
void VecN<N, T>::random(void) 
{
  std::transform(this->begin(), this->end(), this->begin(), Random<T>());
}


template<size_t N, class T>
VecN<N, T> &VecN<N, T>::operator=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), Setval<T>(f));
  return (*this);
}

template<size_t N, class T>
VecN<N, T> &VecN<N, T>::operator+=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        std::bind(std::plus<T>(), std::placeholders::_1, f));
  return (*this);
}

template<size_t N, class T>
VecN<N, T> &VecN<N, T>::operator+=(const VecN<N, T> &w) 
{
  std::transform(this->begin(), this->end(), 
             w.begin(), 
             this->begin(), 
             std::plus<T>() );
  return (*this);
}

template<size_t N, class T>
VecN<N, T> &VecN<N, T>::operator-=(const VecN<N, T> &w)
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::minus<T>() );
  return (*this);
}

template<size_t N, class T>
VecN<N, T> &VecN<N, T>::operator-=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind(std::minus<T>(), std::placeholders::_1, f));
  return (*this);
}

template<size_t N, class T>
VecN<N, T> &VecN<N, T>::operator*=(T f) 
{ 
  std::transform(this->begin(), this->end(), this->begin(), 
        bind(std::multiplies<T>(), std::placeholders::_1, f));
  return *this;  
}

template<size_t N, class T>
VecN<N, T> &VecN<N, T>::operator*=(const VecN<N, T> &w)
{  
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::multiplies<T>() );
  return (*this);  
}

template<size_t N, class T>
VecN<N, T> &VecN<N, T>::operator/=(T f) 
{
  std::transform(this->begin(), this->end(), this->begin(), 
        bind(std::divides<T>(), std::placeholders::_1, f));
  return (*this);
}

template<size_t N, class T>
VecN<N, T> &VecN<N, T>::operator/=(const VecN<N, T> &w) 
{
  std::transform(this->begin(), this->end(), 
        w.begin(), 
        this->begin(), 
        std::divides<T>() );
  return (*this);
}


// Definition of extra functions supporting the VecN class

// find max element

template<size_t N, class T>
T max(const VecN<N, T> &v)
{
  return *std::max_element(v.begin(), v.end());
}

// find min element 

template<size_t N, class T>
T min(const VecN<N, T> &v) 
{
  return *std::min_element(v.begin(), v.end());
}

// the + operator

template<size_t N, class T>
VecN<N, T> operator+(const VecN<N, T> &v, T f) 
{
  VecN<N, T> c(v);
  return c += f;
}

template<size_t N, class T>
VecN<N, T> operator+(T f, const VecN<N, T> &v) 
{
  VecN<N, T> c(v);
  return c += f;
}

template<size_t N, class T>
VecN<N, T> operator+(const VecN<N, T> &v, const VecN<N, T> &w) 
{
  VecN<N, T> c(v);
  return c += w;
}

template<size_t N, class T>
VecN<N, T> operator-(const VecN<N, T> &v, T f)
{
  VecN<N, T> c(v);
  return c -= f;
}

template<size_t N, class T>
VecN<N, T> operator-(T f, const VecN<N, T> &v)
{
  VecN<N, T> c = -v;
  return c += f;
}

template<size_t N, class T>
VecN<N, T> operator-(const VecN<N, T> &v, const VecN<N, T> &w)
{
  VecN<N, T> c(v);
  return c -= w;
}

template<size_t N, class T>
VecN<N, T> operator-(const VecN<N, T> &v)
{
  VecN<N, T> c;
  std::transform(v.begin(), v.end(), c.begin(), std::negate<T>());
  return c;
}

template<size_t N, class T>
T dot(const VecN<N, T> &v, const VecN<N, T> &w) 
{
  return std::inner_product(v.begin(), v.end(), w.begin(), static_cast<T>(0));
}

template<size_t N, class T>
VecN<N, T> operator*(const VecN<N, T> &v, const VecN<N, T> &w) 
{
  VecN<N, T> c(v);
  return c*=w;
}

template<size_t N, class T>
VecN<N, T> operator*(T f, const VecN<N, T> &v)
{
  VecN<N, T> c(f);
  return c *= v;
}

template<size_t N, class T>
VecN<N, T> operator*(const VecN<N, T> &v, T f)
{
  VecN<N, T> c(f);
  return c *= v;
}

template<size_t N, class T>
VecN<N, T> operator/(const VecN<N, T> &v, T f) 
{
  VecN<N, T> c(v);
  return c /= f;
}

template<size_t N, class T>
VecN<N, T> operator/(T f, const VecN<N, T> &v) 
{
  VecN<N, T> c;
  std::transform(v.begin(), v.end(), c.begin(), bind1st(std::divides<T>(), f));
  return c;
}

// element by element division

template<size_t N, class T>
VecN<N, T> operator/(const VecN<N, T> &v, const VecN<N, T> &w) 
{
  VecN<N, T> c(v);
  return c /= w;
}

// Math functions

template<size_t N, class T>
VecN<N, T> sin(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sin ) );
  return b;  
}  

template<size_t N, class T>
VecN<N, T> cos(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( cos ) );
  return b;  
}  

template<size_t N, class T>
VecN<N, T> tan(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( tan ) );
  return b;  
}  

template<size_t N, class T>
VecN<N, T> asin(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( asin ) );
  return b;  
}  

template<size_t N, class T>
VecN<N, T> acos(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( acos ) );
  return b;  
}  

template<size_t N, class T>
VecN<N, T> atan(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( atan ) );
  return b;  
}  

template<size_t N, class T>
VecN<N, T> exp(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( exp ) );
  return b;  
}  

template<size_t N, class T>
VecN<N, T> log(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( log ) );
  return b;  
}  

template<size_t N, class T>
VecN<N, T> sqrt(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( sqrt ) );
  return b;  
}
  
template<size_t N, class T>
VecN<N, T> abs(const VecN<N, T> &v) 
{ 
  VecN<N, T> b;  
  std::transform(v.begin(), v.end(), b.begin(), (T(*)(T))( std::abs) );
  return b;  
}

template<size_t N, class T>
VecN<N, T> pow(const VecN<N, T> &v, T exp) 
{ 
  VecN<N, T> b; 
  std::transform(v.begin(), v.end(), b.begin(), 
                 bind(Power<T>(), std::placeholders::_1, exp));  
  return b;  
}  

template<size_t N, class T>
VecN<N, T> pow(const VecN<N, T> &v, int exp) 
{ 
  // we've got to find a more efficient way to do this...
  return pow(v, static_cast<T>(exp));
}  

// max of 2 vectors

template<size_t N, class T>
VecN<N, T> max(const VecN<N, T> &v, const VecN<N, T> &w) 
{
  VecN<N, T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) > w(i) ? v(i) : w(i);
  return res;
}

template<size_t N, class T> 
VecN<N, T> max(const VecN<N, T> &v, T f) 
{
  VecN<N, T> res;
  std::transform(v.begin(), v.end(), res.begin(), 
                 bind(Max<T>(), std::placeholders::_1, f));
  return res;
}

template<size_t N, class T>
VecN<N, T> max(T f, const VecN<N, T> &v) 
{
  return max(v, f);
}

// min of 2 arguments

template<size_t N, class T> 
VecN<N, T>  min(const VecN<N, T> &v, const VecN<N, T> &w) 
{
  VecN<N, T> res;
  for (size_t i = 0; i < v.size(); ++i)  
    res[i] = v(i) < w(i) ? v(i) : w(i);
  return res;
}

template<size_t N, class T> VecN<N, T>  min(const VecN<N, T> &v, T f) 
{
  VecN<N, T> res;
  std::transform(v.begin(), v.end(), res.begin(), 
                 bind(Min<T>(), std::placeholders::_1, f));
  return res;
}

template<size_t N, class T>
VecN<N, T> min(T f, const VecN<N, T> &v) 
{
  return min(v, f);
}

/*
sum of elements
*/

template<size_t N, class T>
T sum(const VecN<N, T> &v) 
{
  return std::accumulate(v.begin(), v.end(), static_cast<T>(0));
}


// return index of min(v)

template<size_t N, class T>
size_t index_min(VecN<N, T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = min_element(beg, end);
  return static_cast<size_t>(i - beg);
}

// return index of max(v)

template<size_t N, class T>
size_t index_max(VecN<N, T> &v)
{
  T *beg = &v[0];
  T *end = &v[v.size()];  
  T *i = max_element(beg, end);
  return static_cast<size_t>(i - beg);
}

//*******************************************

template <size_t N, class T>
std::ostream& operator<<(std::ostream& s, const VecN<N, T>& v)
{
  for_each(v.begin(), v.end(), Print<T>(s));
  return s;
}

//
// template instantiations
//

template class VecN<2, double>;
template class VecN<3, double>;
template class VecN<4, double>;
template class VecN<6, double>;
template class VecN<9, double>;

template VecN<2, double> operator-(double, const VecN<2, double>&);
template VecN<3, double> operator-(double, const VecN<3, double>&);
template VecN<4, double> operator-(double, const VecN<4, double>&);
template VecN<6, double> operator-(double, const VecN<6, double>&);
template VecN<9, double> operator-(double, const VecN<9, double>&);


template double max(const VecN<2, double>&);
template double max(const VecN<3, double>&);
template double max(const VecN<4, double>&);
template double max(const VecN<6, double>&);
template double max(const VecN<9, double>&);

template double min(const VecN<2, double>&);
template double min(const VecN<3, double>&);
template double min(const VecN<4, double>&);
template double min(const VecN<6, double>&);
template double min(const VecN<9, double>&);

template double dot(const VecN<2, double>&, const VecN<2, double>&);
template double dot(const VecN<3, double>&, const VecN<3, double>&);
template double dot(const VecN<4, double>&, const VecN<4, double>&);
template double dot(const VecN<6, double>&, const VecN<6, double>&);
template double dot(const VecN<9, double>&, const VecN<9, double>&);

template std::ostream& operator<<(std::ostream& s, const VecN<2, double>&);
template std::ostream& operator<<(std::ostream& s, const VecN<3, double>&);
template std::ostream& operator<<(std::ostream& s, const VecN<4, double>&);
template std::ostream& operator<<(std::ostream& s, const VecN<6, double>&);
template std::ostream& operator<<(std::ostream& s, const VecN<9, double>&);

template VecN<2, double> operator+(const VecN<2, double>&, double);
template VecN<3, double> operator+(const VecN<3, double>&, double);
template VecN<4, double> operator+(const VecN<4, double>&, double);
template VecN<6, double> operator+(const VecN<6, double>&, double);
template VecN<9, double> operator+(const VecN<9, double>&, double);

template VecN<2, double> operator+(const VecN<2, double>&, const VecN<2, double>&);
template VecN<3, double> operator+(const VecN<3, double>&, const VecN<3, double>&);
template VecN<4, double> operator+(const VecN<4, double>&, const VecN<4, double>&);
template VecN<6, double> operator+(const VecN<6, double>&, const VecN<6, double>&);
template VecN<9, double> operator+(const VecN<9, double>&, const VecN<9, double>&);

template VecN<2, double> operator-(const VecN<2, double>&, double);
template VecN<3, double> operator-(const VecN<3, double>&, double);
template VecN<4, double> operator-(const VecN<4, double>&, double);
template VecN<6, double> operator-(const VecN<6, double>&, double);
template VecN<9, double> operator-(const VecN<9, double>&, double);

template VecN<2, double> operator-(const VecN<2, double>&, const VecN<2, double>&);
template VecN<3, double> operator-(const VecN<3, double>&, const VecN<3, double>&);
template VecN<4, double> operator-(const VecN<4, double>&, const VecN<4, double>&);
template VecN<6, double> operator-(const VecN<6, double>&, const VecN<6, double>&);
template VecN<9, double> operator-(const VecN<9, double>&, const VecN<9, double>&);

template VecN<2, double> operator*(const VecN<2, double>&, double);
template VecN<3, double> operator*(const VecN<3, double>&, double);
template VecN<4, double> operator*(const VecN<4, double>&, double);
template VecN<6, double> operator*(const VecN<6, double>&, double);
template VecN<9, double> operator*(const VecN<9, double>&, double);

template VecN<2, double> operator*(double, const VecN<2, double>&);
template VecN<3, double> operator*(double, const VecN<3, double>&);
template VecN<4, double> operator*(double, const VecN<4, double>&);
template VecN<6, double> operator*(double, const VecN<6, double>&);
template VecN<9, double> operator*(double, const VecN<9, double>&);

template VecN<2, double> operator*(const VecN<2, double>&, const VecN<2, double>&);
template VecN<3, double> operator*(const VecN<3, double>&, const VecN<3, double>&);
template VecN<4, double> operator*(const VecN<4, double>&, const VecN<4, double>&);
template VecN<6, double> operator*(const VecN<6, double>&, const VecN<6, double>&);
template VecN<9, double> operator*(const VecN<9, double>&, const VecN<9, double>&);

template VecN<2, double> operator/(const VecN<2, double>&, double);
template VecN<3, double> operator/(const VecN<3, double>&, double);
template VecN<4, double> operator/(const VecN<4, double>&, double);
template VecN<6, double> operator/(const VecN<6, double>&, double);
template VecN<9, double> operator/(const VecN<9, double>&, double);

template VecN<2, double> operator/(const VecN<2, double>&, const VecN<2, double>&);
template VecN<3, double> operator/(const VecN<3, double>&, const VecN<3, double>&);
template VecN<4, double> operator/(const VecN<4, double>&, const VecN<4, double>&);
template VecN<6, double> operator/(const VecN<6, double>&, const VecN<6, double>&);
template VecN<9, double> operator/(const VecN<9, double>&, const VecN<9, double>&);

template VecN<2, double> abs(const VecN<2, double>&);
template VecN<3, double> abs(const VecN<3, double>&);
template VecN<4, double> abs(const VecN<4, double>&);
template VecN<6, double> abs(const VecN<6, double>&);
template VecN<9, double> abs(const VecN<9, double>&);

template VecN<2, double> pow(const VecN<2, double>&, double);
template VecN<3, double> pow(const VecN<3, double>&, double);
template VecN<4, double> pow(const VecN<4, double>&, double);
template VecN<6, double> pow(const VecN<6, double>&, double);
template VecN<9, double> pow(const VecN<9, double>&, double);

template VecN<2, double> pow(const VecN<2, double>&, int);
template VecN<3, double> pow(const VecN<3, double>&, int);
template VecN<4, double> pow(const VecN<4, double>&, int);
template VecN<6, double> pow(const VecN<6, double>&, int);
template VecN<9, double> pow(const VecN<9, double>&, int);

template VecN<2, double> exp(const VecN<2, double>&);
template VecN<3, double> exp(const VecN<3, double>&);
template VecN<4, double> exp(const VecN<4, double>&);
template VecN<6, double> exp(const VecN<6, double>&);
template VecN<9, double> exp(const VecN<9, double>&);

template VecN<2, double> log(const VecN<2, double>&);
template VecN<3, double> log(const VecN<3, double>&);
template VecN<4, double> log(const VecN<4, double>&);
template VecN<6, double> log(const VecN<6, double>&);
template VecN<9, double> log(const VecN<9, double>&);

template VecN<2, double> sqrt(const VecN<2, double>&);
template VecN<3, double> sqrt(const VecN<3, double>&);
template VecN<4, double> sqrt(const VecN<4, double>&);
template VecN<6, double> sqrt(const VecN<6, double>&);
template VecN<9, double> sqrt(const VecN<9, double>&);

template VecN<2, double> sin(const VecN<2, double>&);
template VecN<3, double> sin(const VecN<3, double>&);
template VecN<4, double> sin(const VecN<4, double>&);
template VecN<6, double> sin(const VecN<6, double>&);
template VecN<9, double> sin(const VecN<9, double>&);

template VecN<2, double> cos(const VecN<2, double>&);
template VecN<3, double> cos(const VecN<3, double>&);
template VecN<4, double> cos(const VecN<4, double>&);
template VecN<6, double> cos(const VecN<6, double>&);
template VecN<9, double> cos(const VecN<9, double>&);

template VecN<2, double> tan(const VecN<2, double>&);
template VecN<3, double> tan(const VecN<3, double>&);
template VecN<4, double> tan(const VecN<4, double>&);
template VecN<6, double> tan(const VecN<6, double>&);
template VecN<9, double> tan(const VecN<9, double>&);

template VecN<2, double> asin(const VecN<2, double>&);
template VecN<3, double> asin(const VecN<3, double>&);
template VecN<4, double> asin(const VecN<4, double>&);
template VecN<6, double> asin(const VecN<6, double>&);
template VecN<9, double> asin(const VecN<9, double>&);

template VecN<2, double> acos(const VecN<2, double>&);
template VecN<3, double> acos(const VecN<3, double>&);
template VecN<4, double> acos(const VecN<4, double>&);
template VecN<6, double> acos(const VecN<6, double>&);
template VecN<9, double> acos(const VecN<9, double>&);

template VecN<2, double> atan(const VecN<2, double>&);
template VecN<3, double> atan(const VecN<3, double>&);
template VecN<4, double> atan(const VecN<4, double>&);
template VecN<6, double> atan(const VecN<6, double>&);
template VecN<9, double> atan(const VecN<9, double>&);


