#include <vector>

#ifndef MNT_MULTI_ARRAY_ITER
#define MNT_MULTI_ARRAY_ITER

/**
 * A class that stores netcdf attributes
 */

struct MultiArrayIter_t {

  // dimensions along each axis
  std::vector<std::size_t> dims;

  std::vector<std::size_t> prodDims;

  std::size_t bigIndex;

  std::size_t ntot;

};

/**
 * Constructor
 * @param ndims number of dimensions
 * @param dims dimensions
 * @return error code (0 is OK)
 */
extern "C"
int mnt_multiarrayiter_new(MultiArrayIter_t** self, int ndims, const std::size_t* dims);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_multiarrayiter_del(MultiArrayIter_t** self);


/**
 * Set the iterator to the beginning
 * @return error code (0 is OK)
 */
extern "C"
int mnt_multiarrayiter_begin(MultiArrayIter_t** self);

/**
 * Increment the iterator
 * @return error code (0 is OK)
 */
extern "C"
int mnt_multiarrayiter_next(MultiArrayIter_t** self);

/**
 * Get the number of iterations
 * @param n number (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_multiarrayiter_getNumIters(MultiArrayIter_t** self, std::size_t* n);


/**
 * Get the indices
 * @param inds indices (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_multiarrayiter_getIndices(MultiArrayIter_t** self, std::size_t indices[]);

#endif // MULTI_ARRAY_ITER
