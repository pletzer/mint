#include <mntRegridAxis.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#undef NDEBUG // turn on asserts

void testLinear() {
    std::vector<double> x{0., 1., 2., 3., 4., 5. , 6., 7., 8., 9., 10.};
    RegridAxis_t* interp = NULL;
    int ier = mnt_regridaxis_new(&interp);
    assert(ier == 0);

    ier = mnt_regridaxis_build(&interp, (int) x.size(), &x[0]);
    assert(ier == 0);

    double target = 1.1;
    int indices[2];
    double weights[2];
    ier = mnt_regridaxis_getPointWeights(&interp, target, indices, weights);
    assert(ier == 0);
    std::cout << "target: " << target 
              << " interval: " << indices[0] << ',' << indices[1]
              << " weights: " << weights[0] << ',' << weights[1] << '\n';

    ier = mnt_regridaxis_del(&interp);
    assert(ier == 0);
}

void testQuadratic(int n) {
    std::vector<double> t(n);
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        t[i] = double(i);
        x[i] = 1.0 + t[i]*t[i];
    }


    RegridAxis_t* interp = NULL;
    int ier = mnt_regridaxis_new(&interp);
    assert(ier == 0);

    ier = mnt_regridaxis_build(&interp, (int) x.size(), &x[0]);
    assert(ier == 0);

    double target = 1.1;
    int indices[2];
    double weights[2];
    ier = mnt_regridaxis_getPointWeights(&interp, target, indices, weights);
    assert(ier == 0);
    std::cout << "target: " << target 
              << " interval: " << indices[0] << ',' << indices[1]
              << " weights: " << weights[0] << ',' << weights[1] << '\n';

    ier = mnt_regridaxis_del(&interp);
    assert(ier == 0);
}

void testDecreasing(int n) {
    std::vector<double> t(n);
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        t[i] = double(i);
        // decreasing
        x[i] = 1.0 - t[i]*t[i];
    }


    RegridAxis_t* interp = NULL;
    int ier = mnt_regridaxis_new(&interp);
    assert(ier == 0);

    ier = mnt_regridaxis_build(&interp, (int) x.size(), &x[0]);
    assert(ier == 0);

    double target = 0.99;
    int indices[2];
    double weights[2];
    ier = mnt_regridaxis_getPointWeights(&interp, target, indices, weights);
    assert(ier == 0);
    std::cout << "target: " << target 
              << " interval: " << indices[0] << ',' << indices[1]
              << " weights: " << weights[0] << ',' << weights[1] << '\n';

    ier = mnt_regridaxis_del(&interp);
    assert(ier == 0);
}

void testCells(int n, const double xtargets[]) {
    std::vector<double> t(n);
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        t[i] = double(i);
        // decreasing
        x[i] = 1.0 + t[i];
    }


    RegridAxis_t* interp = NULL;
    int ier = mnt_regridaxis_new(&interp);
    assert(ier == 0);

    ier = mnt_regridaxis_build(&interp, (int) x.size(), &x[0]);
    assert(ier == 0);

    double indexBounds[2];
    ier = mnt_regridaxis_getCellIndexBounds(&interp, xtargets, indexBounds);
    assert(ier == 0);

    // conservative interpolation weights
    std::cout << "target interval: " << xtargets[0] << ',' << xtargets[1] 
              << " index bounds = " << indexBounds[0] << ',' << indexBounds[1] << '\n';
    for (int i = std::floor(indexBounds[0]); i < std::floor(indexBounds[1]) + 1; ++i) {
        std::cout << " cell index = " << i << " weight = " 
                  << std::min(indexBounds[1], double(i + 1)) - std::max(indexBounds[0], double(i)) << '\n';
    }

    ier = mnt_regridaxis_del(&interp);
    assert(ier == 0);
}


int main(int argc, char** argv) {

    testLinear();
    testQuadratic(10);
    testDecreasing(10);

    // target cell is inside
    double xtargets[] = {1.1, 1.2};
    testCells(10, xtargets);

    // target cell matches a grid cell
    xtargets[0] = 1.; xtargets[1] = 2.;
    testCells(10, xtargets);

    // target cell overlaps several grid cells
    xtargets[0] = 1.1; xtargets[1] = 5.9;
    testCells(10, xtargets);

    std::cout << "Success\n";

    return 0;
}
