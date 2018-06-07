
/**
 * @brief Parse command line arguments
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */

#include "CmdLineArgParser.h"

template <class T>
T CmdLineArgParser::get(const std::string& name) const {
  std::cerr << "ERROR. Invalid type, use get<TYPE>(...)\n";
}

template <>
double CmdLineArgParser::get(const std::string& name) const {
  double res = -std::numeric_limits<double>::max();
  std::map<std::string, double>::const_iterator i = this->doubleArg.find(name);
  if (i != this->doubleArg.end()) {
    res = i->second;
  }
  return res;
}

template <>
int CmdLineArgParser::get(const std::string& name) const {
  int res = -std::numeric_limits<int>::max();
  std::map<std::string, int>::const_iterator i = this->intArg.find(name);
  if (i != this->intArg.end()) {
    res = i->second;
  }
  return res;
}
  
template <>
std::string CmdLineArgParser::get(const std::string& name) const {
  std::string res = "";
  std::map<std::string, std::string>::const_iterator
    i = this->stringArg.find(name);
  if (i != this->stringArg.end()) {
    res = i->second;
  }
return res;
}
  

template <>
bool CmdLineArgParser::get(const std::string& name) const {
  bool res = false;
  std::map<std::string, bool>::const_iterator
    i = this->boolArg.find(name);
  if (i != this->boolArg.end()) {
    res = i->second;
  }
  return res;
}
