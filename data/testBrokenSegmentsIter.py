from math import pi
from cubedsphere_reader import CubedsphereReader
from broken_line_iter import BrokenLineIter
from broken_segments_iter import BrokenSegmentsIter

def test(filename, points):

    csr = CubedsphereReader(filename=filename)
    bl = BrokenLineIter(points)
    bs = BrokenSegmentsIter(csr.getUnstructuredGrid(), 
                            csr.getUnstructuredGridCellLocator(), 
                            bl)
    count = 0
    for s in bs:
        cellId = s.getCellId()
        xia = s.getBegCellParamCoord()
        xib = s.getEndCellParamCoord()
        ta = s.getBegLineParamCoord()
        tb = s.getEndLineParamCoord()
        print('seg {} in cell {} t = {} -> {} xi = {} -> {}'.format(count, cellId, ta, tb, xia, xib))
        count += 1

    totalT = s.getIntegratedParamCoord()
    print('Integrated t = {}'.format(totalT))
    assert(abs(totalT - 1.0) < 1.e-10)


if __name__ == '__main__':
    print('Test1')
    test(filename='mesh_C4.nc', points=[(0., 1.1),(1., 1.1)])
    print('Test2')
    test(filename='mesh_C4.nc', points=[(0., 0.1),(2., 0.1)])