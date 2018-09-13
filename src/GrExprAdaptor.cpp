// -*-c++-*-
// $Id: GrExprAdaptor.cpp 524 2013-10-01 03:24:44Z pletzer $

// grin includes
#include <GrExprAdaptor.h>

// standard includes
#include <deque>
#include <set>

GrExprAdaptor::GrExprAdaptor(const std::string& expr) {

  // squeeze all blanks out 
  for (size_t i = 0; i < expr.size(); ++i) {
    char c = expr[i];
    if (c != ' ') {
      this->expr.push_back(c);
    }
  }

  // define some default operators
  this->defineUnaryOperator("sin");
  this->defineUnaryOperator("cos");
  this->defineUnaryOperator("tan");
  this->defineUnaryOperator("log");
  this->defineUnaryOperator("exp");
  this->defineUnaryOperator("-");

  this->defineBinaryOperator("+", 0);
  this->defineBinaryOperator("*", 1);
  //this->defineBinaryOperator("-", 0);
  this->defineBinaryOperator("/", 1);

  // break the string into tokens
  this->tokenize();
}

GrExprAdaptor::~GrExprAdaptor() {
}

void
GrExprAdaptor::defineUnaryOperator(const std::string& name) {
  this->unaryOps.insert(name);
  this->tokenize();
}

void
GrExprAdaptor::defineBinaryOperator(const std::string& name, 
                                   int priority) {
  this->binaryOps.insert( std::pair<std::string, int>(name, priority) );
  this->tokenize();
}

std::string
GrExprAdaptor::getPrefixExpr() const {

  std::deque<std::string> opStack;
  std::deque<std::string> exprList;

  // number of arguments, either 2 for binary operation, 1 for unary and 0 if element 
  // is a variable
  std::deque<int> numArgList;

  for (std::deque<std::string>::const_reverse_iterator 
       it = this->tokens.rbegin(); 
       it != this->tokens.rend(); ++it) {
    const std::string& elem = *it;
    if (this->binaryOps.find(elem) != this->binaryOps.end()) {
      // it's a binary operator. check if there already is another
      // operator in opStack. 
      if (opStack.size() > 0) {
        // move otherOp to exprList
        // add elem to opStack
        std::string otherOp = opStack[0];
        std::map<std::string, int>::const_iterator 
                              itBinOtherOp = this->binaryOps.find(otherOp);
        std::map<std::string, int>::const_iterator 
                              itBinElemOp = this->binaryOps.find(elem);
        if (itBinOtherOp != this->binaryOps.end() &&
            itBinElemOp != this->binaryOps.end() &&
            itBinOtherOp->second >= itBinElemOp->second) {
          opStack.pop_front();
          if (otherOp != ")") {
            exprList.push_front(otherOp);
            numArgList.push_front(2); // binary op
          }
        }
      }
      opStack.push_front(elem);
    }
    else if (elem == ")") {
      opStack.push_front(elem);
    }
    else if (elem == "(") {
      // pour all the operators into exprList until ")" is encountered, 
      // in reverse order
      bool foundRightParen = false;
      while (opStack.size() > 0 && !foundRightParen) {
        std::string otherOp = opStack[0];
        opStack.pop_front();
        if (otherOp != ")") {
          exprList.push_front(otherOp);
          numArgList.push_front(2); // binary op
        }
        else {
          foundRightParen = true;
        }
      }
    }
    else if (elem != "") {
      // variable or unary operator
      exprList.push_front(elem);
      if (this->unaryOps.find(elem) != this->unaryOps.end()) {
        numArgList.push_front(1); // unary op
      }
      else {
        numArgList.push_front(0); // a variable has no argument
      }
    }
  }

  // push the remaining operators in the stack
  while (opStack.size() > 0) {
    std::string otherOp = opStack[0];
    opStack.pop_front();
    exprList.push_front(otherOp);
    numArgList.push_front(2); // binary op
  }

  // push the + operator to the front if the first element of exprList
  // is a variable. 
  if (numArgList.size() > 0 && numArgList[0] == 0) {
    exprList.push_front("+");
    if (numArgList.size() > 1 && numArgList[1] == 0) {
      numArgList.push_front(2);
    }
    else {
      // unary op
      numArgList.push_front(1);
    }
  }

  // reduce the expression list by identifying all the binary and
  // unary operators and condensing the arguments into a single 
  // string
  bool hasReduced = true;
  while (hasReduced) {
    hasReduced = this->reduce(exprList, numArgList);
  }
  
  // convert the list into a string
  std::string res;
  for (size_t i = 0; i < exprList.size(); ++i) {
    
    // token in the expression list
    std::string elem = exprList[i];

    // space delimiter
    if (res.size() > 0) {
      res += " ";
    }

    // add element
    res += elem;
  }

  // replace all commas by a space
  for (size_t i = 0; i < res.size(); ++i) {
    if (res[i] == ',') {
      res[i] = ' ';
    }
  }

  return res;
}

bool
GrExprAdaptor::reduce(std::deque<std::string>& exprList, 
                      std::deque<int>& numArgList) const {

  bool hasReduced = false;
  if (exprList.size() <= 1) {
    return hasReduced;
  }

  for (int i = (int)exprList.size() - 1; i >= 0; --i) {

    if (numArgList[i] == 2) {
      // binary operator
      exprList[i] = "(" + exprList[i] + ' ' 
        + exprList[i + 1] + ' ' + exprList[i + 2] + ")";
      exprList.erase(exprList.begin() + i + 1, 
                     exprList.begin() + i + 3);
      numArgList[i] = 0;
      numArgList.erase(numArgList.begin() + i + 1,
                       numArgList.begin() + i + 3);
      hasReduced = true;
      break;
    }

    if (numArgList[i] == 1) {
      // unary operator
      exprList[i] = "(" + exprList[i] + ' ' 
        + exprList[i + 1] + ")";
      exprList.erase(exprList.begin() + i + 1, 
                     exprList.begin() + i + 2);
      numArgList[i] = 0;
      numArgList.erase(numArgList.begin() + i + 1,
                       numArgList.begin() + i + 2);
      hasReduced = true;
      break;
    }
  }
  
  return hasReduced;
}

void
GrExprAdaptor::tokenize() {

  this->tokens.clear();
  
  // list of delimiters for each atomic object, be it 
  // a variable, an operator, or a parenthesis

  // left delimiters
  std::set<std::string> leftDelimiters;
  leftDelimiters.insert("(");
  leftDelimiters.insert(")");
  // add the unary operators to the list of left delimiters
  for (std::set<std::string>::const_iterator 
         it = this->unaryOps.begin(); 
       it != this->unaryOps.end(); ++it) {
    leftDelimiters.insert(*it);
  }
  // add the binary operators to the list of left delimiters
  for (std::map<std::string, int>::const_iterator 
         it = this->binaryOps.begin(); 
       it != this->binaryOps.end(); ++it) {
    leftDelimiters.insert(it->first);
  }
  
  // right delimiters
  std::set<std::string> rightDelimiters;
  rightDelimiters.insert(")");
  // add the binary operators to the list of right delimiters
  for (std::map<std::string, int>::const_iterator 
         it = this->binaryOps.begin(); 
       it != this->binaryOps.end(); ++it) {
    rightDelimiters.insert(it->first);
  }
  
  size_t lpos = 0;
  size_t rpos = this->expr.size();
  while (rpos > 0 && lpos < rpos) {
 
    // find the right most left delimiter
    std::string lDelim("");
    for (std::set<std::string>::const_iterator lit = leftDelimiters.begin();
         lit != leftDelimiters.end(); ++lit) {
      const std::string& d = *lit;
      size_t pos = this->expr.rfind(d, rpos - 1);
      if (pos != std::string::npos) {
        if (pos >= lpos) {
          lpos = pos;
          lDelim = d;
        }
      }
    }
    
    // start position of the token
    size_t pos1 = lpos + lDelim.size();
    size_t n = rpos - pos1;

    // get the token
    std::string token = this->expr.substr(pos1, n);

    // look for right delimiters in the token. Further break the
    // token into smaller token until no more right delimiters 
    // are found.
    size_t rp = token.size();
    bool foundRightDelimiter = true;
    while (foundRightDelimiter) {
      foundRightDelimiter = false;
      for (std::set<std::string>::const_iterator rit = rightDelimiters.begin();
           rit != rightDelimiters.end(); ++rit) {
        const std::string& d = *rit;
        size_t pos = token.rfind(d, rp - 1);
        if (pos != std::string::npos) {
          // found a right delimiter
          this->tokens.push_front( token.substr(pos, rp) );
          rp = pos;
          foundRightDelimiter = true;
          break;
        }
      }
    }
   
    // insert the remaining token into the list
    this->tokens.push_front( token.substr(0, rp) );

    // reset the left/right positions
    lpos = 0;
    rpos = pos1;
    
    // push the delimiter
    if (lDelim != "" && lDelim != " ") {
      this->tokens.push_front(lDelim);
      rpos -= lDelim.size();
    }
  }
  
}
