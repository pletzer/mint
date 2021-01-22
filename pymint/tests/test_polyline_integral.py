from mint import Grid
from mint import PolylineIntegral
import numpy
import sys
from pathlib import Path


def test_simple():

    # create the grid
    gr = Grid()

    # 2 cells
    points = numpy.array([(0., 0., 0.),
                          (1., 0., 0.),
                          (1., 1., 0.),
                          (0., 1., 0.),
                          (1., 0., 0.),
                          (2., 0., 0.),
                          (2., 1., 0.),
                          (1., 1., 0.)]).reshape(2, 4, 3)
    gr.setPoints(points)

    pli = PolylineIntegral()
    pli.setGrid(gr)
    xyz = numpy.array([(0., 0., 0.),
    	               (2., 0., 0.),
    	               (2., 1., 0.),
    	               (0., 1., 0.)])
    pli.setPolyline(xyz)

    data = numpy.array([(1., 2., 1., 2.), 
    	                (1., 2., 1., 2.)])

    flux = pli.getIntegral(data)
    print(f'total flux: {flux:.3f}')
    assert abs(flux - 2.0) < 1.e-10


if __name__ == '__main__':

    data_dir = Path('./data')
    if len(sys.argv) >= 2:
        data_dir = Path(sys.argv[1])

    test_simple()


