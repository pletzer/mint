from mint import MultiArrayIter
import numpy
from pathlib import Path
from tempfile import TemporaryDirectory

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def test_create_destroy():
    mai = MultiArrayIter([])
    mai.begin()
    num_iter = mai.getNumIters()
    assert(num_iter == 0)

def test_1():
    mai = MultiArrayIter([1])
    mai.begin()
    num_iter = mai.getNumIters()
    assert(num_iter == 1)

def test_10():
    mai = MultiArrayIter([10])
    mai.begin()
    num_iter = mai.getNumIters()
    assert(num_iter == 10)
    count = 0
    for i in range(num_iter):
        inds = mai.getIndices()
        print(f'{i} indices: {inds}')
        mai.next()
        count += 1
    assert(count == num_iter)

def test_3d():
    dims = numpy.array([2, 3, 4])
    mai = MultiArrayIter(dims=dims)
    mai.begin()
    num_iter = mai.getNumIters()
    assert(num_iter == numpy.prod(dims))
    count = 0
    for i in range(num_iter):
        inds = mai.getIndices()
        assert( numpy.all([inds[i] < dims[i] for i in range(len(dims))]) )
        print(f'{i} indices: {inds}')
        mai.next()
        count += 1
    assert(count == num_iter)


