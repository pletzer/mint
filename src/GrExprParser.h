// -*-c++-*-
// $Id: GrExprParser.h 519 2013-09-30 00:03:02Z pletzer $

#ifndef GR_EXPR_PARSER_H
#define GR_EXPR_PARSER_H

// standard includes
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <cmath>

#include <MvVector.h>
typedef Vector<double> Vec;

/**
 * @brief Simple expression parser for Polish prefix notation, e.g.
 * (+ a b) for a + b
 */

class GrExprParser {

  typedef Vec* (*vectFunc_t)(const std::vector<Vec*>&);
  typedef Vec* (*unaryFunc_t)(const Vec*, void*);

public:

  /**
   * Constructor
   * @param n vector argument length
   * @param expr expression in prefix notation, e.g. (+ a b)
   * @note spaces are delimitors, the first string is the 
   * operator
   */
  GrExprParser(size_t n, const std::string& expr);
  
  /**
   * Destructor
   */
  virtual ~GrExprParser();

  /**
   * Define a constant
   * @param name name of the constant
   * @param val value
   */
  void defineConstant(const std::string& name, double val);

  /**
   * Define a variable
   * @param name name of the object in the expression
   * @param object pointer to the object, object is owned by caller
   */
  void defineVariable(const std::string& name, Vec* object);

  /**
   * Define an operator
   * @param name name of operator, as it appears in expression
   * @param fct function
   */
  void defineOperator(const std::string& name, vectFunc_t fct);
  void defineOperator(const std::string& name, unaryFunc_t fct, 
                      void* params = 0);

  /**
   * Evaluate expression
   * @return pointer to the result, this class owns the output object
   * @note will return null pointer if an error was encountered
   */
  Vec* eval();
  
  /** 
   * Convert expresion from infix to prefix notation
   * @param inExpr expression in infix notation, e.g. 3*sin(x+y*z)
   * @return expression in prefix notation, e.g. (* 3 (sin (+ x (* y z))))
   */
  std::string convertInfixToPrefixNotation(const std::string& inExpr) const;

 private:
  
  /**
   * Fill in list of variable names appearing in expression
   * @param varnames on input an empty list, will contain the list of 
   *                   variables on output
   * @param opnames on input an empty list, will contain the list of 
   *                   operators on output
   * @param expr on input an empty string, will contain the reduced expression
   *               on output
   * @note this method should be called before eval.
   */
  void fillInVarNamesFromExpr(std::set<std::string>& varnames, 
                              std::set<std::string>& opnames,
                              std::string& expr) const;
  
  /**
   * Remove any duplicate space in an expression
   * @param expr expression
   * @return expression
   */
  std::string squeeze(const std::string& expr) const;

  /**
   * Tokenize a string using space as delimiter
   * @param expr expression
   * @return array of tokens
   */
  std::vector<std::string> tokenize(const std::string& expr) const;

  /** 
   * Addition operator
   * @param args arguments
   * @return args[0] + args[1] + args[2]...
   */
  static Vec* add(const std::vector<Vec*>& args) {
    Vec* res = 0;
    if (args.size() > 0) {
      size_t n = args[0]->size(); // assume all other args have the same size
      res = new Vec(n, 0.0);
      for (size_t i = 0; i < args.size(); ++i) {
        *res += *args[i];
      }
    }
    return res;
  }

  /** 
   * Subtraction operator
   * @param args arguments
   * @return args[0] - args[1] - args[2]...
   */
  static Vec* sub(const std::vector<Vec*>& args) {
    Vec* res = 0;
    if (args.size() > 0) {
      // assume all other args have the same size
      res = new Vec(*(args[0]));
      for (size_t i = 1; i < args.size(); ++i) {
        *res -= *args[i];
      }
    }
    return res;
  }

  /**
   * Multiplication operator
   * @param args arguments
   * @return args[0] * args[1] * args[2]...
   */
  static Vec* mul(const std::vector<Vec*>& args) {
    Vec* res = 0;
    if (args.size() > 0) {
      size_t n = args[0]->size(); // assume all other args have the same size
      res = new Vec(n, 1.0);
      for (size_t i = 0; i < args.size(); ++i) {
         *res *= *args[i];
      }
    }
    return res;
  }

  /**
   * Division operator
   * @param args arguments
   * @return args[0] / args[1] / args[2]...
   */
  static Vec* div(const std::vector<Vec*>& args) {
    Vec* res = 0;
    if (args.size() > 0) {
      // assume all other args have the same size
      res = new Vec(*(args[0]));
      for (size_t i = 1; i < args.size(); ++i) {
         *res /= *args[i];
      }
    }
    return res;
  }

  /**
   * Negative (unary)
   * @param x argument
   * @param params (not used)
   * @return -x elementwise
   */
  static Vec* neg(const Vec* x, void* params = 0) {
    Vec* res = new Vec( -(*x) );
    return res;
  }
  
  /**
   * Inverse (unary)
   * @param x argument
   * @param params (not used)
   * @return 1/x elementwise
   */
  static Vec* inv(const Vec* x, void* params = 0) {
    Vec* res = new Vec( 1.0/(*x) );
    return res;
  }
  
  /**
   * Square root (unary)
   * @param x argument
   * @param params (not used)
   * @return sqrt elementwise
   */
  static Vec* sqrt(const Vec* x, void* params = 0) {
    size_t n = x->size();
    Vec* res = new Vec(n);
    for (size_t i = 0; i < n; ++i) {
      (*res)(i) = std::sqrt( (*x)(i) );
    }
    return res;
  }
  
  /**
   * Natural log (unary)
   * @param x argument
   * @param params (not used)
   * @return log(x) elementwise
   */
  static Vec* log(const Vec* x, void* params = 0) {
    size_t n = x->size();
    Vec* res = new Vec(n);
    for (size_t i = 0; i < n; ++i) {
      (*res)(i) = std::log( (*x)(i) );
    }
    return res;
  }
  
  /**
   * Log in base 10 (unary)
   * @param x argument
   * @param params (not used)
   * @return log10(x) elementwise
   */
  static Vec* log10(const Vec* x, void* params = 0) {
    size_t n = x->size();
    Vec* res = new Vec(n);
    for (size_t i = 0; i < n; ++i) {
      (*res)(i) = std::log10( (*x)(i) );
    }
    return res;
  }

  /**
   * sin (unary)
   * @param x argument
   * @param params (not used)
   * @return sin(x) elementwise
   */
  static Vec* sin(const Vec* x, void* params = 0) {
    size_t n = x->size();
    Vec* res = new Vec(n);
    for (size_t i = 0; i < n; ++i) {
      (*res)(i) = std::sin( (*x)(i) );
    }
    return res;
  }

  /**
   * cos (unary)
   * @param x argument
   * @param params (not used)
   * @return sin(x) elementwise
   */
  static Vec* cos(const Vec* x, void* params = 0) {
    size_t n = x->size();
    Vec* res = new Vec(n);
    for (size_t i = 0; i < n; ++i) {
      (*res)(i) = std::cos( (*x)(i) );
    }
    return res;
  }

  /**
   * tan (unary)
   * @param x argument
   * @param params (not used)
   * @return sin(x) elementwise
   */
  static Vec* tan(const Vec* x, void* params = 0) {
    size_t n = x->size();
    Vec* res = new Vec(n);
    for (size_t i = 0; i < n; ++i) {
      (*res)(i) = std::tan( (*x)(i) );
    }
    return res;
  }

  /** Holds pointers to the inputs */
  std::map<std::string, Vec*> input;

  /** Holds pointers to the temporaries */
  std::map<std::string, Vec*> temps;

  /** Expression after applying reduction operations */
  std::string expr;

  /** map of function pointers */
  std::map<std::string, vectFunc_t> ops;

  /** map of unary function pointers */
  std::map<std::string, unaryFunc_t> unOps;

  /** map of parameters corresponding to each of the above unary fct */
  std::map<std::string, void*> params;

  /** constants */
  std::map<std::string, double> consts;

  /** list of variable names */
  std::vector<std::string> varnames;

  /** size of vector arguments */
  size_t n;
};

#endif // GR_EXPR_PARSER_H
