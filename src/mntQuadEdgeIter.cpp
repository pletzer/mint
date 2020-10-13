#include <mntQuadEdgeIter.h>

QuadEdgeIter::QuadEdgeIter() {

    // create a VTK QUAD cell and extract the edge to node connectivity
    this->cell = vtkQuad::New(); // 2d
    for (size_t i = 0; i < cell->GetNumberOfPoints(); ++i) {
        this->cell->GetPointIds()->SetId(i, i);
    }

    this->allParamCoords = this->cell->GetParametricCoords();

    // get the edges of the cell and correct for the direction - we want the 
    // edges always to point in the positive direction
    for (int edgeId = 0; edgeId < cell->GetNumberOfEdges(); ++edgeId) {

        const vtkIdType* i01 = this->cell->GetEdgeArray(edgeId);
        int i0 = i01[0];
        int i1 = i01[1];

        // check if the direction is positive
        double dir = 0;
        double* xi0 = &(this->allParamCoords)[3*i0];
        double* xi1 = &(this->allParamCoords)[3*i1];
        for (size_t d = 0; d < 3; ++d) {
            // dir is either -1 (negative direction), 0 or +1 (positive direction)
            dir += xi1[d] - xi0[d];
        }
        if (dir < 0) {
            // swap the direction of the edge so that the edge always points to the 
            // positive direction
            int i0Bis = i0;
            i0 = i1;
            i1 = i0Bis;
        }

        this->edgeIds.push_back(edgeId);
        this->iBegs.push_back(i0);
        this->iEnds.push_back(i1);
    }

}

QuadEdgeIter::~QuadEdgeIter() {
    this->cell->Delete();
}

int 
QuadEdgeIter::getNumberOfEdges() const {
    return (int) this->edgeIds.size();
}

void
QuadEdgeIter::getCellPointIds(int edgeId, int* iBeg, int* iEnd) const {
    *iBeg = this->iBegs[edgeId];
    *iEnd = this->iEnds[edgeId];
}

void
QuadEdgeIter::getParamCoords(int edgeId, double** xiBeg, double** xiEnd) {
    int i0 = this->iBegs[edgeId];
    int i1 = this->iEnds[edgeId];
    *xiBeg = &(this->allParamCoords)[3*i0];
    *xiEnd = &(this->allParamCoords)[3*i1];
}
