// -*-c++-*-
// $Id: GrExprParser.cpp 495 2013-09-24 21:01:45Z pletzer $

// grin includes
#include <GrExprParser.h>

// standard includes
#include <sstream>
#include <cstdlib> // strtod
#include <deque>

GrExprParser::GrExprParser(size_t n, const std::string& expr) {
  this->n = n;
  this->expr = this->squeeze(expr);

  // some default constants
  this->defineConstant("pi", 3.141592653589793116);
  this->defineConstant("e", 2.718281828459045091);

  // set multi-ops
  this->defineOperator("+", GrExprParser::add);
  this->defineOperator("-", GrExprParser::sub);
  this->defineOperator("*", GrExprParser::mul);
  this->defineOperator("/", GrExprParser::div);

  // set unary ops
  this->defineOperator("sqrt", GrExprParser::sqrt);
  this->defineOperator("log", GrExprParser::log);
  this->defineOperator("log10", GrExprParser::log10);
  this->defineOperator("sin", GrExprParser::sin);
  this->defineOperator("cos", GrExprParser::cos);
  this->defineOperator("tan", GrExprParser::tan);
  this->defineOperator("-", GrExprParser::neg);
}

GrExprParser::~GrExprParser() {
  std::map<std::string, Vec*>::iterator it;
  for (it = this->temps.begin(); it != this->temps.end(); ++it) {
    if (it->second) {
      delete it->second;
      it->second = 0;
    }
  }
}

void
GrExprParser::defineConstant(const std::string& name, double val) {
  // create a "temporary" 
  Vec* temp = new Vec(this->n, val);
  this->temps.insert( std::pair<std::string, Vec*>(name, temp) );
}

void
GrExprParser::defineVariable(const std::string& name, Vec* object) {
  this->input.insert( std::pair<std::string, Vec*>(name, object) ) ;
}

Vec* 
GrExprParser::eval() {

  size_t lpos = this->expr.rfind('(');
  size_t rpos = this->expr.find(')', lpos);
  if (lpos == std::string::npos || 
      rpos == std::string::npos) {
    std::ostringstream error;
    error << "GrExprParser::eval ERROR: unbalanced parentheses in ";
    error << this->expr << '\n';
    throw error.str();
  }
  std::string subExpr = this->expr.substr(lpos+1, rpos-lpos-1);
  std::vector<std::string> tokens = this->tokenize(subExpr);

  // fill in the operator argument list
  std::vector<Vec*> args;
  // first element is operator
  for (size_t i = 1; i < tokens.size(); ++i) {
    std::map<std::string, Vec*>::const_iterator it;
    // look into the input var container
    it = this->input.find(tokens[i]);
    if (it != this->input.end()) {
      args.push_back(it->second);
    }
    else {
      // is this a temporary object?
      it = this->temps.find(tokens[i]);
      if (it != this->temps.end()) {
        args.push_back(it->second);
      }
      else {
        // are we dealing with a numeral?
        char* endPtr;
        double val = strtod(tokens[i].c_str(), &endPtr);
        if (!*endPtr) {
          // turn the number into a Vec, store in the list
          // of temporaries as this will allow the instance
          // to clean up upon destruction. The name of the 
	        // variable is the string version of the constant.
          Vec* temp = new Vec(this->n, 1);
          *temp *= val;
          this->temps.insert( std::pair<std::string, Vec*>(tokens[i], temp) );
          args.push_back(temp);
        }
	else {
          std::ostringstream error;
	  error << "GrExprParser::eval ERROR: unknown object named ";
	  error << tokens[i] << " in sub expr " << subExpr << '\n';
	  throw error.str();
	}
      }
    } 
  }

  // create temporary, these are named $0, $1, ...
  Vec* temp = 0;
  std::ostringstream oss;
  oss << "$" << this->temps.size();
  std::string tempName = oss.str();

  const std::string& op = tokens[0];

  if (args.size() == 1) {
    // unary operator
    std::map<std::string, unaryFunc_t>::const_iterator it;
    it = this->unOps.find(op);
    if (it != this->unOps.end()) {
      // matched a unary operator
      temp = it->second(args[0], this->params[op]);
    }
  }

  if (!temp) {
    // multi-var operator
    std::map<std::string, vectFunc_t>::const_iterator it;
    it = this->ops.find(op);
    if (it != this->ops.end()) {
      // matched a multi-var operator
      temp = it->second(args);
    }
  }

  if (temp) {
    // everything's fine...
    std::pair<std::string, Vec*> p(tempName, temp);
    this->temps.insert(p);
  }
  else {
    std::ostringstream error;
    error << "GrExprParser::eval ERROR: undefined operator ";
    error << op << '\n';
    throw error.str();
  }

  // reduce the expression
  this->expr = this->expr.substr(0, lpos-0) + 
    tempName + this->expr.substr(rpos+1);

  // check if this is the last expression to evaluate
  if (this->expr.rfind('(') == std::string::npos) {
    return temp;
  }
  else {
    // otherwise recursively call this method
    return this->eval();
  }
}

std::string
GrExprParser::squeeze(const std::string& expr) const {
  std::string res;
  bool skip = true;
  for (size_t i = 0; i < expr.size(); ++i) {
    char c = expr[i];
    if (c != ' ') {
      res.push_back(c);
      if (c != '(' && c != ')') {
        skip = false;
      }
    }
    else if (!skip) {
      res.push_back(c);
      skip = true;
    }
  }
  return res;
}

std::vector<std::string> 
GrExprParser::tokenize(const std::string& array) const {
  // "+ ab cdef" -> {"+", "ab", "cdef"}
  std::vector<std::string> res;
  std::string elem = "";
  for (size_t i = 0; i < array.size(); ++i) {
    char c = array[i];
    if (c != ' ') {
      elem.push_back(c);
    }
    if (c == ' ' || i == array.size() - 1) {
      if (elem != "") {
        res.push_back(elem);
      }
      elem = "";
    }
  }
  return res;
}

void
GrExprParser::defineOperator(const std::string& name, 
                                   vectFunc_t fct) {
  std::pair<std::string, vectFunc_t> p(name, fct);
  this->ops.insert(p);
}

void
GrExprParser::defineOperator(const std::string& name, 
                                   unaryFunc_t fct,
                                   void* params) {
  std::pair<std::string, unaryFunc_t> p(name, fct);
  this->unOps.insert(p);
  this->params.insert(std::pair<std::string, void*>(name, params));
}

void
GrExprParser::fillInVarNamesFromExpr(std::set<std::string>& varnames,
                                             std::set<std::string>& opnames,
                                             std::string& expr) const {

  // this method should be called BEFORE eval;
  if (expr == "") {
    expr = this->expr;
  }

  size_t lpos = expr.rfind('(');
  size_t rpos = expr.find(')', lpos);
  if (lpos == std::string::npos || 
      rpos == std::string::npos) {
    std::ostringstream error;
    error << "GrExprParser::fillInVarNamesFromExpr ERROR: ";
    error << "unbalanced parentheses in " << expr << '\n';
    throw error.str();
  }
  std::string subExpr = expr.substr(lpos+1, rpos-lpos-1);
  std::vector<std::string> tokens = this->tokenize(subExpr);

  // is the operator already defiend?
  if (opnames.find(tokens[0]) == opnames.end()) {
    opnames.insert(tokens[0]);
  }

  // first element is operator
  for (size_t i = 1; i < tokens.size(); ++i) {
    if (varnames.find(tokens[i]) == varnames.end()) {
      // only add to list new variable names
      varnames.insert(tokens[i]);
    }
  }

  // reduce the expression, cutting out the part we
  // already looked at
  expr = expr.substr(0, lpos-0) + expr.substr(rpos+1);

  // check if this is the last expression to evaluate
  if (expr.rfind('(') == std::string::npos) {
    // done
    return;
  }
  else {
    // otherwise recursively call this method
    this->fillInVarNamesFromExpr(varnames, opnames, expr);
  }

}

std::string
GrExprParser::convertInfixToPrefixNotation(const std::string& inExpr) const {
  
  // Iterate over each element from the back. An element is always bounded
  // by: start/end of the string, any of the supported operator (+, -, *,
  // ...), and space. An element can either be a variable or an operator. 
  // Operators are either binary (a + b) or unary (sin(x)). 
  // Variable goes to the back of exprList, operators go into the front of
  // opStack. When a new operator is added to opStack which already 
  // contains an operator, the former operator goes into exprList if 
  // this operator has precedence over the earlier stored operator, 
  // otherwise it is added to the opStack. At the end, all the operators
  // are taken out of the opStack in the inverse order they were put 
  // in. 

  std::vector<std::string> unaryOps;
  unaryOps.push_back("sin");
  unaryOps.push_back("cos");
  unaryOps.push_back("-");

  // binary operators and their priority
  std::map<std::string, int> binaryOps;
  binaryOps.insert( std::pair<std::string, int>("+", 0) );
  binaryOps.insert( std::pair<std::string, int>("-", 0) );
  binaryOps.insert( std::pair<std::string, int>("*", 1) );
  binaryOps.insert( std::pair<std::string, int>("/", 1) );

  std::vector<std::string> varDelimiters;
  varDelimiters.push_back(" ");
  varDelimiters.push_back("(");
  varDelimiters.push_back(")");
  for (size_t i = 0; i < unaryOps.size(); ++i) {
    varDelimiters.push_back(unaryOps[i]);
  }
  for (std::map<std::string, int>::const_iterator 
    it = binaryOps.begin(); it != binaryOps.end(); ++it) {
    varDelimiters.push_back(it->first);
  }

  // tokenize the expression in reverse order
  std::vector<std::string> tokens;
  size_t begPos = 0;
  size_t endPos = inExpr.size() - 1;
  while (begPos <= endPos) {
    // find the position of the last delimiter
    std::string delimiter = "";
    for (size_t j = 0; j < varDelimiters.size(); ++j) {
      size_t bp = inExpr.rfind(varDelimiters[j], endPos);
      if (bp != std::string::npos && bp > begPos) {
        begPos = bp;
        delimiter = varDelimiters[j];
      }
    }
    // extract the variable or operator
    std::string varOrOp = inExpr.substr(begPos + delimiter.size(), endPos);
    // remove any spaces
    std::string varOrOpNoSpaces = "";
    for (size_t j = 0; j < varOrOp.size(); ++j) {
      char c = varOrOp[j];
      if (c != ' ') {
        varOrOpNoSpaces.push_back(c);
      }
    }
    // now add to list
    tokens.push_back(varOrOpNoSpaces);
    endPos = begPos;
    begPos = 0;
    if (delimiter != " ") {
      tokens.push_back(delimiter);
    }
  }

  // create the prefix expression, note that tokens is
  // already in reverse order
  std::deque<std::string> opStack;
  std::deque<std::string> exprList;
  for (size_t i = 0; i < tokens.size(); ++i) {
    std::string token = tokens[i];
    // determine if token is a variable, a delimiter, or an operator
    // ( and ) are moved to the opStack
    // unary and variables are pushed to exprList
    bool isOpenParen = false;
    bool isCloseParen = false;
    bool isBinaryOp = false;
    if (token == ")") {
      isCloseParen = true;
    }
    else if (token == "(") {
      isOpenParen = true;
    }
    else {
      for (std::map<std::string, int>::const_iterator 
             it = binaryOps.begin(); it != binaryOps.end(); ++it) {
        if (token == it->first) {
          isBinaryOp = true;
        }
      }
    }
    // if not a binary op, a closing or opening parenthesis, then 
    // it must be a unary or a variable
    bool isVarOrUnaryOp = (!isOpenParen) && (!isCloseParen) && (!isBinaryOp);
    if (isVarOrUnaryOp) {
      exprList.push_front(token);
    }
    else if (isCloseParen) {
      opStack.push_front(token);
    }
    else if (isBinaryOp) {
      // if there already is an operator in the stack 
      // and if this operator is stronger, then pull 
      // this operator out, then push the new operator.
      if (opStack.size() > 0 && binaryOps[opStack[0]] > binaryOps[token]) {
        // e.g., token = "+" and opStack[0] = "*"
        exprList.push_front(opStack[0]);
        opStack.pop_front();
      }
      opStack.push_front(token);
    }
    else {
      // opening parenthesis
      // pour the operators in inverse order they came in
      // until ")" is found
      bool rightParenFound = false;
      while (opStack.size() > 0 && !rightParenFound) {
        if (opStack[0] == ")") {
          rightParenFound = true;
        }
        else {
          exprList.push_front(opStack[0]);
        }
        opStack.pop_front();
      }
    }
  }

  // Build the prefix expression. Put brackets around unary and
  // binary operations
  std::string preExpr;
  for (size_t i = 0; i < exprList.size(); ++i) {
    preExpr += exprList[i];
  }

  return preExpr;
}
