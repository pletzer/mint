// -*-C++-*- 
// $Id: MvFunctors.h 180 2013-06-14 12:54:31Z pletzer $

#ifndef _FUNCTORS_H_
#define _FUNCTORS_H_

#include <functional>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <cstdlib>

// List of functors for Vector and Matrix classes

template<class T> struct Random : public std::unary_function<T, T> 
{
  Random() : max(static_cast<T>(RAND_MAX)) {}
  T operator() (T x) {return static_cast<T>(rand())/max;}
  T max;
};

template<class T> struct Setval : public std::unary_function<T, T>
{
  Setval(T fin) : f(fin) {};
  T operator() (T x) {return f;}
  T f;
};

template<class T> struct Power : public std::binary_function<T, T, T> 
{
  T operator() (T x, T y) const { return std::pow(x, y); }
};

template<class T> struct Max : public std::binary_function<T, T, T> 
{
  T operator() (T x, T y) const { return std::max<T>(x, y); }
};

template<class T> struct Min : public std::binary_function<T, T, T> 
{
  T operator() (T x, T y) const { return std::min<T>(x, y); }
};

template<class T> struct Print : public std::unary_function<T, void>
{
  std::ostream& os;
  Print(std::ostream& out) : os(out) {
    os.precision(13);
  }
  void operator() (T x) { os << x << ' '; }
};

#endif
