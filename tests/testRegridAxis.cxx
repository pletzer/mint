#include <mntRegridAxis.h>
#include <vector>
#include <cmath>
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

void testCellInside(int n) {
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

    double xtargets[] = {1.1, 1.2};
    int numCellWeights;
    ier = mnt_regridaxis_getNumCellWeights(&interp, xtargets, &numCellWeights);
    assert(ier == 0);
    assert(numCellWeights == 1);

    int indices[numCellWeights];
    double weights[numCellWeights];
    ier = mnt_regridaxis_getCellWeights(&interp, indices, weights);
    assert(ier == 0);
    std::cout << "target interval: " << xtargets[0] << ',' << xtargets[1] << " indices/weights: ";
    for (int i = 0; i < numCellWeights; ++i) {
        std::cout << indices[i] << '/' << weights[0];
    }
    std::cout << '\n';
    //assert(indices[0] == 1);
    //assert(std::fabs(weights[0] - 0.1) < 1.e-10);

    ier = mnt_regridaxis_del(&interp);
    assert(ier == 0);
}


int main(int argc, char** argv) {

    testLinear();
    testQuadratic(10);
    testDecreasing(10);
    testCellInside(10);

    std::cout << "Success\n";

    return 0;
}
