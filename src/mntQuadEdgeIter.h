#include <vector>
#include <vtkQuad.h>

#ifndef MNT_QUAD_EDGE_ITER
#define MNT_QUAD_EDGE_ITER

class QuadEdgeIter {

public:

    /**
     * Constructor
     */
    QuadEdgeIter();

    /**
     * Destructor
     */
    ~QuadEdgeIter();

    /**
     * Get the number of edges
     * @return number
     */
    int getNumberOfEdges() const; 

    /**
     * Get the two cell point Ids for given edge
     * @param edgeId edge index
     * @param iBeg index of begin point in range [0...num cell points( (output)
     * @param iEnd index of end point in range [0...num cell points( (output)
     */
    void getCellPointIds(int edgeId, int* iBeg, int* iEnd) const;

    /**
     * Get the two cell point Ids for given edge
     * @param edgeId edge index
     * @param xiBeg start point in parametric space (output)
     * @param xiEnd iend point in parametric space (output)
     */
    void getParamCoords(int edgeId, double** xiBeg, double** xiEnd);

private:

    // maps an edge index to a pair of vertex indices
    vtkQuad* cell;
    double* allParamCoords;
    std::vector<int> edgeIds;
    std::vector<int> iBegs;
    std::vector<int> iEnds;
};

#endif // MNT_QUAD_EDGE_ITER
