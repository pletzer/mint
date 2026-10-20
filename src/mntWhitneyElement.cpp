#include <mntWhitneyElement.h>
#include <mntWeights.h>
#include <mntGlobal.h>


void
WhitneyElement::getCovariantComponents(double xi, double eta, const double dofs[],
                                       double* a, double* b) const {
    double as[MNT_MAX_DOFS_PER_CELL];
    double bs[MNT_MAX_DOFS_PER_CELL];
    this->evalBasis(xi, eta, as, bs);
    *a = 0.0;
    *b = 0.0;
    for (int i = 0; i < this->getNumDofs(); ++i) {
        *a += dofs[i]*as[i];
        *b += dofs[i]*bs[i];
    }
}


void
WhitneyElementOrder1::evalBasis(double xi, double eta, double a[], double b[]) const {
    // south
    a[0] = 1.0 - eta;
    b[0] = 0.0;
    // east
    a[1] = 0.0;
    b[1] = xi;
    // north
    a[2] = eta;
    b[2] = 0.0;
    // west
    a[3] = 0.0;
    b[3] = 1.0 - xi;
}

void
WhitneyElementOrder1::getCovariantComponents(double xi, double eta, const double dofs[],
                                             double* a, double* b) const {
    // closed form of the generic sum, skips the zero terms
    *a = dofs[0]*(1.0 - eta) + dofs[2]*eta;
    *b = dofs[3]*(1.0 - xi) + dofs[1]*xi;
}

void
WhitneyElementOrder1::integrateBasis(const Vec3& xia, const Vec3& xib, double w[]) const {
    for (int i = 0; i < MNT_NUM_EDGES_PER_QUAD; ++i) {
        const double* xi0 = &QUAD_POSITIVEXI_EDGE_DIRECTION[6*i];
        const double* xi1 = &QUAD_POSITIVEXI_EDGE_DIRECTION[6*i + 3];
        w[i] = computeWeight(xi0, xi1, xia, xib);
    }
}

int
WhitneyElementOrder1::getCounterclockwiseSign(int dofIndex) const {
    // the north and west edges point clockwise in the positive xi/eta orientation
    return (dofIndex < 2)? 1: -1;
}

int
WhitneyElementOrder1::getCellDofs(Grid_t* grid, vtkIdType cellId, const double data[], int placement,
                                  double dofs[]) const {
    int ier = 0;
    if (placement == MNT_CELL_BY_CELL_DATA) {
        for (int i = 0; i < MNT_NUM_EDGES_PER_QUAD; ++i) {
            dofs[i] = data[cellId*MNT_NUM_EDGES_PER_QUAD + i];
        }
    }
    else {
        std::size_t edgeId;
        int edgeSign;
        for (int i = 0; i < MNT_NUM_EDGES_PER_QUAD; ++i) {
            if (mnt_grid_getEdgeId(&grid, cellId, i, &edgeId, &edgeSign) != 0) {
                dofs[i] = 0.0;
                ier++;
                continue;
            }
            dofs[i] = data[edgeId] * edgeSign;
        }
    }
    return ier;
}
