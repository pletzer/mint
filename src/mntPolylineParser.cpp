#include "mntPolylineParser.h"
#include <cstdlib>

Vector<double> 
PolylineParser::parsePosition(const std::string& posStr) const {
    Vector<double> res(this->ndims, 0);
    size_t posBeg = 0;
    size_t posEnd = posStr.find(',');
    size_t counter = 0;
    while (posEnd != std::string::npos && counter < this->ndims) {
        double number = std::atof( posStr.substr(posBeg, posEnd - posBeg - 1).c_str() );
        res[counter] = number;
        counter++;
    }
    return res;
}

void
PolylineParser::parse(const std::string& expr) {
    this->points.clear();
    const char leftDelim = '(';
    const char rghtDelim = ')';
    size_t leftPos = 0;
    size_t rghtPos = expr.size();
    while (true) {
        leftPos = expr.find(leftDelim, leftPos);
        rghtPos = expr.find(rghtDelim, leftPos);
        if (leftPos == std::string::npos || rghtPos == std::string::npos) {
            // could not find a point or parentheses don't match
            break;
        }
        leftPos++;
        size_t n = rghtPos - leftPos;
        std::string posStr = expr.substr(leftPos, n);
        this->points.push_back(this->parsePosition(posStr));
    }
}

