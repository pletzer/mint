#include <mntWhitneyElement.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>

// edge end points in parameter space, positive xi/eta direction
const double EDGES[] = {
    0., 0., 0.,   1., 0., 0., // south
    1., 0., 0.,   1., 1., 0., // east
    0., 1., 0.,   1., 1., 0., // north
    0., 0., 0.,   0., 1., 0.  // west
};

Vec3 point(const double* p) {
    Vec3 res;
    for (std::size_t d = 0; d < 3; ++d) res[d] = p[d];
    return res;
}

/**
 * Integral of each basis function along xia -> xib using the midpoint rule
 */
void quadratureIntegrals(const WhitneyElement& elem, const Vec3& xia, const Vec3& xib,
                         double w[]) {
    const int n = 1000;
    const int numDofs = elem.getNumDofs();
    double a[MNT_MAX_DOFS_PER_CELL], b[MNT_MAX_DOFS_PER_CELL];
    for (int i = 0; i < numDofs; ++i) w[i] = 0.0;
    for (int k = 0; k < n; ++k) {
        double t = (k + 0.5)/n;
        double xi = xia[0] + t*(xib[0] - xia[0]);
        double eta = xia[1] + t*(xib[1] - xia[1]);
        elem.evalBasis(xi, eta, a, b);
        for (int i = 0; i < numDofs; ++i) {
            w[i] += (a[i]*(xib[0] - xia[0]) + b[i]*(xib[1] - xia[1]))/n;
        }
    }
}

void testDuality(const WhitneyElement& elem) {
    // the integral of basis function i along edge j is delta_ij
    double w[MNT_MAX_DOFS_PER_CELL];
    for (int j = 0; j < elem.getNumDofs(); ++j) {
        elem.integrateBasis(point(&EDGES[6*j]), point(&EDGES[6*j + 3]), w);
        for (int i = 0; i < elem.getNumDofs(); ++i) {
            double expected = (i == j)? 1.0: 0.0;
            assert(std::abs(w[i] - expected) < 1.e-14);
        }
        // reversing the edges that run clockwise gives the counterclockwise sign
        elem.integrateBasis(point(&EDGES[6*j + 3]), point(&EDGES[6*j]), w);
        assert(std::abs(w[j] + 1.0) < 1.e-14);
    }
    std::cout << "testDuality: OK\n";
}

void testIntegrateBasis(const WhitneyElement& elem) {
    // exact segment integrals agree with quadrature of the basis functions
    const double segs[][6] = {
        {0.1, 0.2, 0.,   0.9, 0.7, 0.},
        {0.8, 0.1, 0.,   0.3, 0.95, 0.},
        {0.5, 0.5, 0.,   0.5, 0.0, 0.},
        {0.0, 0.3, 0.,   1.0, 0.3, 0.},
    };
    double w[MNT_MAX_DOFS_PER_CELL], wq[MNT_MAX_DOFS_PER_CELL];
    for (const auto& s : segs) {
        elem.integrateBasis(point(&s[0]), point(&s[3]), w);
        quadratureIntegrals(elem, point(&s[0]), point(&s[3]), wq);
        for (int i = 0; i < elem.getNumDofs(); ++i) {
            assert(std::abs(w[i] - wq[i]) < 1.e-6);
        }
    }
    std::cout << "testIntegrateBasis: OK\n";
}

void testCovariantComponents(const WhitneyElement& elem) {
    // the specialized getCovariantComponents matches the generic sum over basis functions
    const double dofs[] = {1.3, -0.7, 2.1, 0.4, 0.9, -1.1, 0.2, 0.5, -0.3, 1.7, 0.8, -2.0};
    for (double xi = 0.; xi <= 1.0; xi += 0.125) {
        for (double eta = 0.; eta <= 1.0; eta += 0.125) {
            double a, b, aGeneric, bGeneric;
            elem.getCovariantComponents(xi, eta, dofs, &a, &b);
            elem.WhitneyElement::getCovariantComponents(xi, eta, dofs, &aGeneric, &bGeneric);
            assert(std::abs(a - aGeneric) < 1.e-14);
            assert(std::abs(b - bGeneric) < 1.e-14);
        }
    }
    std::cout << "testCovariantComponents: OK\n";
}

void testCounterclockwiseSigns(const WhitneyElement& elem) {
    const int expected[] = {1, 1, -1, -1};
    for (int i = 0; i < elem.getNumDofs(); ++i) {
        assert(elem.getCounterclockwiseSign(i) == expected[i]);
    }
    std::cout << "testCounterclockwiseSigns: OK\n";
}

int main(int argc, char** argv) {

    WhitneyElementOrder1 elem;
    assert(elem.getOrder() == 1);
    assert(elem.getNumDofs() == 4);

    testDuality(elem);
    testIntegrateBasis(elem);
    testCovariantComponents(elem);
    testCounterclockwiseSigns(elem);

    return 0;
}
