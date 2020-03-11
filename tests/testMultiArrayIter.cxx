#include <netcdf.h>
#include <string>
#include <iostream>
#undef NDEBUG // turn on asserts
#include <cassert>
#include "mntMultiArrayIter.h"

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
            std::cout << ' ' << inds[i] << ' ';
        }
        std::cout << '\n';
        ier = mnt_multiarrayiter_next(&mai);
        assert(ier == 0);
    }

    ier = mnt_multiarrayiter_del(&mai);
    assert(ier == 0);
}

int main(int argc, char** argv) {

    std::vector<size_t> dims0;
	test(dims0);

    return 0;
}
