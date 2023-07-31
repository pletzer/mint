import numpy as np
from pathlib import Path
import mint

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

def test_1():
    nf = mint.NcFieldRead(DATA_DIR / 'cs8_wind.nc', 'u_in_w2h')
    assert nf.getNumDims() == 1
    nedges = nf.getDim(iAxis=0)

    dimName = nf.getDimName(iAxis=0)
    assert nf.getDimName(iAxis=0) == 'ncs_edge'

    data = nf.data()
    assert abs(data.sum() - 603.5819343269815 < 1.e-10)

if __name__ == '__main__':
    test_1()