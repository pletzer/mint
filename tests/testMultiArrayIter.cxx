#include <netcdf.h>
#include <string>
#include <iostream>
#undef NDEBUG // turn on asserts
#include <cassert>
#include "mntMultiArrayIter.h"

void test0() {

    // tests the case ndims == 0
    // expect one iteration
    // getIndices should be a no operation

    int ier;

    MultiArrayIter_t* mai;

    ier = mnt_multiarrayiter_new(&mai, 0, NULL);
    assert(ier == 0);

    size_t niters;
    ier = mnt_multiarrayiter_getNumIters(&mai, &niters);
    assert(ier == 0);
    assert(niters == 1);

    ier = mnt_multiarrayiter_del(&mai);
    assert(ier == 0);
}

void test(const std::vector<size_t>& dims) {

    int ier;

    MultiArrayIter_t* mai;
    
    int ndims = (int) dims.size();
    ier = mnt_multiarrayiter_new(&mai, ndims, &dims[0]);
    assert(ier == 0);

    size_t niters;
    ier = mnt_multiarrayiter_getNumIters(&mai, &niters);
    assert(ier == 0);


    std::vector<size_t> inds(ndims);
    for (size_t i = 0; i < niters; ++i) {

        ier = mnt_multiarrayiter_getIndices(&mai, &inds[0]);
        assert(ier == 0);

        std::cout << i << ": ";
        for (int j = 0; j < ndims; ++j) {
            std::cout << ' ' << inds[j] << ' ';
        }
        std::cout << '\n';

        ier = mnt_multiarrayiter_next(&mai);
        assert(ier == 0);
    }

    ier = mnt_multiarrayiter_del(&mai);
    assert(ier == 0);
}

int main(int argc, char** argv) {

    std::vector<size_t> dims0{};
    test0();

    std::vector<size_t> dims1{2};
    test(dims1);

    std::vector<size_t> dims2{2, 3};
    test(dims2);

    std::vector<size_t> dims3{2, 3, 4};
    test(dims3);

    return 0;
}
