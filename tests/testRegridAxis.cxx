#include <mntRegridAxis.h>
#include <vector>
#include <cassert>
#undef NDEBUG // turn on asserts

void testSimple() {
    std::vector<double> x{1., 5., 6.};
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

int main(int argc, char** argv) {

    testSimple();

    std::cout << "Success\n";

    return 0;
}
