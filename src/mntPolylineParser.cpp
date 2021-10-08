#include "mntLogger.h"
#include "mntPolylineParser.h"
#include <cstdlib>
#include <iostream>

Vec3
PolylineParser::parsePosition(const std::string& posStr) const {
    Vec3 res(0.0);
    std::size_t posBeg, posEnd;
    posBeg = 0;
    for (std::size_t i = 0; i < this->ndims; ++i) {
        // index of ',' or end of string, starting at posBeg
        posEnd = posStr.find(',', posBeg);
        // convert substring to number
        double number = std::atof( posStr.substr(posBeg, posEnd - posBeg).c_str() );
        res[i] = number;
        if (posEnd == std::string::npos && i < (this->ndims) - 1) {
            std::string msg = "missing comma in \"" + posStr + "\"";
            mntlog::warn(__FILE__, __func__, __LINE__, msg);
            break;
        }
        posBeg = posEnd + 1;
    }
    return res;
}

void
PolylineParser::parse(const std::string& expr) {
    this->points.clear();
    const char leftDelim = '(';
    const char rghtDelim = ')';
    std::size_t leftPos = 0;
    std::size_t rghtPos = expr.size();
    while (true) {
        leftPos = expr.find(leftDelim, leftPos);
        rghtPos = expr.find(rghtDelim, leftPos);
        if (leftPos == std::string::npos || rghtPos == std::string::npos) {
            // could not find a point or parentheses don't match
            break;
        }
        leftPos++;
        std::size_t n = rghtPos - leftPos;
        std::string posStr = expr.substr(leftPos, n);
        this->points.push_back(this->parsePosition(posStr));
    }
}

void
PolylineParser::print() const {
    std::cout << "Points:\n";
    for (std::size_t i = 0; i < this->points.size(); ++i) {
        const Vec3& point = this->points[i];
        for (std::size_t j = 0; j < this->ndims; ++j) {
            std::cout << point[j] << ',';
        }
        std::cout << '\n';
    }

}

