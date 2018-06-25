// -*-c++-*-
// $Id: GrExprAdaptor.h 519 2013-09-30 00:03:02Z pletzer $

#ifndef GR_EXPR_ADAPTOR_H
#define GR_EXPR_ADAPTOR_H

// standard includes
#include <string>
#include <set>
#include <deque>
#include <map>
#include <iostream>

/**
 * @brief A class that takes an expression in infix notation
 * and produces a new expression in prefix notation
 */
class GrExprAdaptor {

public:

  /**
   * Constructor
   * @param expr expression in infix notation, e.g. a + b
   */
  GrExprAdaptor(const std::string& expr);
  
  /**
   * Destructor
   */
  virtual ~GrExprAdaptor();

  /**
   * Define a unary operator
   * @param name name of operator
   */
  void defineUnaryOperator(const std::string& name);

  /**
   * Define a binary operator
   * @param name name of operator
   * @param priority priority, large values mean operator takes precedence
   *                 over other operators
   */
  void defineBinaryOperator(const std::string& name, int priority = 0);

  /**
   * Define a multivariable operator
   * @param name name of operator
   * @param numArgs number of arguments
   */
  void defineMultiVarOperator(const std::string& name, int numArgs);

  /**
   * Get the equivalent prefix expression
   */
  std::string getPrefixExpr() const;

 private:
  
  /**
   * Tokenize the expression
   */
  void tokenize();

  /**
   * Reduce an expression
   * @param exprList list of atomic objects (either operators or variables).
   *        On output: a shorted list
   * @param numArgList list of number of argument (1 or 2) for operators, 
   *        0 for variables.
   */
  bool reduce(std::deque<std::string>& exprList, 
              std::deque<int>& numArgList) const;

  /** original expression */
  std::string expr;

  /** unary operators */
  std::set<std::string> unaryOps;

  /** binary operators */
  std::map<std::string, int> binaryOps;

  /** list of tokens */
  std::deque<std::string> tokens;
};

#endif // GR_EXPR_ADAPTOR_H
