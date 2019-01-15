
/**
 * @brief Parse command line arguments
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */

#include "CmdLineArgParser.h"

void 
CmdLineArgParser::set(const std::string& name, double defaultVal, const std::string& help) {
    this->doubleArg[name] = defaultVal;
    this->doubleArgHelp[name] = help;
}
  
void 
CmdLineArgParser::set(const std::string& name, int defaultVal, const std::string& help) {
    this->intArg[name] = defaultVal;
    this->intArgHelp[name] = help;
}
  
void 
CmdLineArgParser::set(const std::string& name, const std::string& defaultVal, const std::string& help) {
    this->stringArg[name] = defaultVal;
    this->stringArgHelp[name] = help;
}
  
void 
CmdLineArgParser::set(const std::string& name, bool defaultVal, const std::string& help) {
    this->boolArg[name] = defaultVal;
    this->boolArgHelp[name] = help;    
}

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
  
bool 
CmdLineArgParser::parse(int argc, char *argv[]) {

  bool isArgVal;
  std::string nm;
  this->execName = argv[0];

  // check if all the options supplied are valid
  for (int a = 1; a < argc; ++a) {
    std::string arg(argv[a]);
    // is arg the name of an option or a value?
    // assume that options are of the form "-X"
    bool isOptionName = (arg[0] == '-')
      && (arg.size() >= 2) && !isdigit(arg[1]);
    if (isOptionName) {
      bool isValidOptName =
        doubleArg.find(arg) != doubleArg.end()
        || intArg.find(arg) != intArg.end()
        || stringArg.find(arg) != stringArg.end()
        || boolArg.find(arg) != boolArg.end();
      if (!isValidOptName) {
        std::cerr << "Warning: '" << arg << "' is not a valid option.\n";
        return false;
      }
    }
  }

  // doubles
  for (std::map<std::string, double>::iterator
       i = this->doubleArg.begin();
       i != this->doubleArg.end(); ++i) {
    std::string optKey = i->first;

    isArgVal = false;
    nm = "";
    for (int a = 1; a < argc; ++a) {
      std::string arg(argv[a]);
      if (isArgVal) {
        i->second = atof(arg.c_str());
        isArgVal = false;
      }
      if (optKey == arg) {
        isArgVal = true;
        nm = optKey;
      }
    }
  }

  // integers
  for (std::map<std::string, int>::iterator
        i = this->intArg.begin();
        i != this->intArg.end(); ++i) {
    std::string optKey = i->first;

    isArgVal = false;
    nm = "";
    for (int a = 1; a < argc; ++a) {
      std::string arg(argv[a]);
      if (isArgVal) {
        i->second = atoi(arg.c_str());
        isArgVal = false;
      }
      if (optKey == arg) {
        isArgVal = true;
        nm = optKey;
      }
    }
  }
    
  // strings
  for (std::map<std::string, std::string>::iterator
       i = this->stringArg.begin();
       i != this->stringArg.end(); ++i) {
    std::string optKey = i->first;

    isArgVal = false;
    nm = "";
    for (int a = 1; a < argc; ++a) {
      std::string arg(argv[a]);
      if (isArgVal) {
        i->second = arg;
        isArgVal = false;
      }
      if (optKey == arg) {
        isArgVal = true;
        nm = optKey;
      }
    }
  }

  // booleans
  for (std::map<std::string, bool>::iterator
       i = this->boolArg.begin();
       i != this->boolArg.end(); ++i) {
    std::string optKey = i->first;
    for (int a = 1; a < argc; ++a) {
      std::string arg(argv[a]);
      if (arg == i->first) {
        // flip the value
        i->second = !i->second;
      }
    }
  }

  return true;
}
  

void 
CmdLineArgParser::help() const {
  std::cout << this->purpose << std::endl;
  std::cout << this->execName << " [options]\n";
  std::cout << "Usage:\n";
  for (std::map<std::string, double>::const_iterator
         i = this->doubleArg.begin();
       i != this->doubleArg.end(); ++i) {
    std::string hlp = this->doubleArgHelp.find(i->first)->second;
    std::cout << "\t" << i->first << " <double#> " << hlp
              << " (" << i->second << ")\n";
  }
  for (std::map<std::string, int>::const_iterator
         i = this->intArg.begin();
       i != this->intArg.end(); ++i) {
    std::string hlp = this->intArgHelp.find(i->first)->second;
    std::cout << "\t" << i->first << " <int#> " << hlp
              << " (" << i->second << ")\n";
  }
  for (std::map<std::string, std::string>::const_iterator
         i = this->stringArg.begin();
       i != this->stringArg.end(); ++i) {
    std::string hlp = this->stringArgHelp.find(i->first)->second;
    std::cout << "\t" << i->first << " <string> " << hlp
              << " (" << i->second << ")\n";
  }
  for (std::map<std::string, bool>::const_iterator
         i = this->boolArg.begin();
       i != this->boolArg.end(); ++i) {
    std::string hlp = this->boolArgHelp.find(i->first)->second;
    std::cout << "\t" << i->first << ' ' << hlp
              << " (" << i->second << ")\n";
  }
  std::cout << this->footer << std::endl;    
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
