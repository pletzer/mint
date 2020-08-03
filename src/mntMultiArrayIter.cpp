#include <mntMultiArrayIter.h>
#include <iostream>

extern "C"
int mnt_multiarrayiter_new(MultiArrayIter_t** self, int ndims, const size_t* dims) {

    *self = new MultiArrayIter_t();
    (*self)->dims.resize(ndims);
    (*self)->prodDims.resize(ndims);
    (*self)->ntot = 1;
    for (int i = 0; i < ndims; ++i) {
        (*self)->dims[i] = dims[i];
        (*self)->ntot *= dims[i];
    }

    // row major
    if (ndims > 0) {
        (*self)->prodDims[0] = 1;
    }
    for (int i = 1; i < ndims; ++i) {
        (*self)->prodDims[i] = (*self)->prodDims[i - 1] * dims[i - 1];
    }

    (*self)->bigIndex = 0;

    return 0;
}

extern "C"
int mnt_multiarrayiter_del(MultiArrayIter_t** self) {
    delete *self;
    return 0;
}


extern "C"
int mnt_multiarrayiter_begin(MultiArrayIter_t** self) {
    (*self)->bigIndex = 0;
    return 0;
}

extern "C"
int mnt_multiarrayiter_next(MultiArrayIter_t** self) {
    (*self)->bigIndex++;
    return 0;
}

extern "C"
int mnt_multiarrayiter_getNumIters(MultiArrayIter_t** self, size_t* n) {
    *n = (*self)->ntot;
    return 0;
}

extern "C"
int mnt_multiarrayiter_getIndices(MultiArrayIter_t** self, size_t indices[]) {
    for (size_t i = 0; i < (*self)->dims.size(); ++i) {
        indices[i] = (*self)->bigIndex / (*self)->prodDims[i] % (*self)->dims[i];
    }
    return 0;
}
