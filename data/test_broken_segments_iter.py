from math import pi
from ugrid_reader import UgridReader
from broken_line_iter import BrokenLineIter
from broken_segments_iter import BrokenSegmentsIter

def test(filename, points):

    ur = UgridReader(filename=filename)
    bl = BrokenLineIter(points)
    bs = BrokenSegmentsIter(ur.getUnstructuredGrid(), 
                            ur.getUnstructuredGridCellLocator(), 
                            bl)
    count = 0
    for s in bs:
        cellId = s.getCellId()
        xia = s.getBegCellParamCoord()
        xib = s.getEndCellParamCoord()
        ta = s.getBegLineParamCoord()
        tb = s.getEndLineParamCoord()
        coeff = s.getCoefficient()
        print('seg {} in cell {} t = {} -> {} xi = {} -> {} coeff = {}'.format(count, cellId, ta, tb, xia, xib, coeff))
        count += 1

    totalT = s.getIntegratedParamCoord()
    print('Integrated t = {}'.format(totalT))
    assert(abs(totalT - 1.0) < 1.e-10)


if __name__ == '__main__':
    print('Test1')
    test(filename='mesh_C4.nc', points=[(0., 1.1),(1., 1.1)])
    print('Test2')
    test(filename='mesh_C4.nc', points=[(0., 0.1),(2., 0.1)])
    print('Test3')
    test(filename='mesh_C4.nc', points=[(0., 1.1),(2., 1.1)])
    print('Test4')
    test(filename='mesh_C4.nc', points=[(0., 1.),(2*pi, 1.)])
    print('Test5')
    test(filename='mesh_C16.nc', points=[(0., 1.4),(2*pi, 1.4)])
    print('Test6')
    test(filename='mesh_C4.nc', points=[(0., 1.),(0., -1.),(2*pi, -1.)])
    print('Test7')
    test(filename='mesh_C4.nc', points=[(3.9269908169872414, -0.6154797086703874,),(3.5342917352885173, -0.7458526607730738,)])
    print('Test8')
    test(filename='mesh_C4.nc', points=[(3.5342917352885173, 0.7458526607730737,),(3.9269908169872414, 0.6154797086703873,)])
    print('Test9')
    test(filename='mesh_C4.nc', points=[(0.7853981633974483, -1.04089353704597, 0.0), (0.3926990816987242, -0.7458526607730738, 0.0)])