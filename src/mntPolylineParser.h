#include <vector>
#include <string>
#include <mntVecN.h>

#ifndef MNT_POLYLINE_PARSER
#define MNT_POLYLINE_PARSER

class PolylineParser {

public:

    /**
     * Constructor
     * @param ndims number of dimensions
     */
    PolylineParser(size_t ndims) {
        this->ndims = ndims;
    }

    /**
     * Parse the expression string
     * @param expr for instance "(1,2,3), (4.5, 6.7, -8.e-9)"
     */
    void parse(const std::string& expr);

    /**
     * Get vector of points
     * @return vector of points
     */
    const std::vector< Vec3 >& getPoints() const {
        return this->points;
    }

    /**
     * Print the points
     */
    void print() const; 


private:

    Vec3 parsePosition(const std::string& posStr) const;

    size_t ndims;
    std::string expr;
    std::vector< Vec3 > points;

};

#endif // MNT_POLYLINE_PARSER
