#include <mntWhitneyElement.h>
#include <mntGlobal.h>


// edge start and end points in parameter space, edges point in the positive xi/eta direction
static const double QUAD_POSITIVEXI_EDGE_DIRECTION[] = {
    0., 0., 0.,   1., 0., 0., // south
    1., 0., 0.,   1., 1., 0., // east
    0., 1., 0.,   1., 1., 0., // north
    0., 0., 0.,   0., 1., 0.  // west
};

/**
 * Integrate the lowest order basis function of an edge along a straight segment
 * @param srcXi0 start point of the edge, in parametric coordinates
 * @param srcXi1 end point of the edge, in parametric coordinates
 * @param xia start point of the segment, in parametric coordinates
 * @param xib end point of the segment, in parametric coordinates
 * @return integral of the edge's basis function along the segment
 */
static double computeEdgeWeight(const double srcXi0[], const double srcXi1[],
                                const Vec3& xia, const Vec3& xib) {

    double weight = 1.0;
    double sgn = 0.0;

    // iterate over the two dimensions of the quad
    for (std::size_t d = 0; d < 2; ++d) { // 2d 

        double xiM = 0.5*(xia[d] + xib[d]);
        double dxi = xib[d] - xia[d];

        // mid point of edge in parameter space
        double xm = 0.5*(srcXi0[d] + srcXi1[d]);
        sgn += srcXi1[d] - srcXi0[d]; // should be either -1, 0 or 1

        // use Lagrange interpolation to evaluate the basis function integral for
        // any of the 3 possible x values in {0, 0.5, 1}. This formula will make 
        // it easier to extend the code to 3d
        double xm00 = xm;
        double xm05 = xm - 0.5;
        double xm10 = xm - 1.0;
        double lag00 = + 2.0 * xm05 * xm10;
        double lag05 = - 4.0 * xm00 * xm10;
        double lag10 = + 2.0 * xm00 * xm05;

        weight *= (1.0 - xiM)*lag00 + dxi*lag05 + xiM*lag10;
    }

    return sgn * weight;
}


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


int
WhitneyElement::getCellDofs(Grid_t* grid, vtkIdType cellId, const double data[], int placement,
                            double dofs[]) const {
    int numFailures = 0;
    std::size_t index;
    int sign;
    for (int i = 0; i < this->getNumDofs(); ++i) {
        if (this->getDofDataIndex(grid, cellId, i, placement, &index, &sign) != 0) {
            dofs[i] = 0.0;
            numFailures++;
            continue;
        }
        dofs[i] = data[index] * sign;
    }
    return numFailures;
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
        w[i] = computeEdgeWeight(xi0, xi1, xia, xib);
    }
}

int
WhitneyElementOrder1::getCounterclockwiseSign(int dofIndex) const {
    // the north and west edges point clockwise in the positive xi/eta orientation
    return (dofIndex < 2)? 1: -1;
}

int
WhitneyElementOrder1::getDofDataIndex(Grid_t* grid, vtkIdType cellId, int dofIndex, int placement,
                                      std::size_t* index, int* sign) const {
    // one degree of freedom per edge
    if (placement == MNT_CELL_BY_CELL_DATA) {
        *index = cellId*MNT_NUM_EDGES_PER_QUAD + dofIndex;
        *sign = 1;
        return 0;
    }
    return mnt_grid_getEdgeId(&grid, cellId, dofIndex, index, sign);
}
