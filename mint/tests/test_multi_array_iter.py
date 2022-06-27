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


