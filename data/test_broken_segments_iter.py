from math import pi
from ugrid_reader import UgridReader
from latlon_reader import LatLonReader
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

def testUm(filename, points, padding):

    ur = LatLonReader(filename=filename, padding=padding)
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

    # UGRID data
    print('Test0 ' + '-'*40)
    eps = 1.73654365e-10
    test(filename='mesh_C4.nc', points=[(pi/2., 0.-eps, 0.),(pi/2. + 2.*pi/16., 0.+eps, 0.)])
    print('Test1 ' + '-'*40)
    test(filename='mesh_C4.nc', points=[(0., 1.1, 0.),(1., 1.1, 0.)])
    print('Test2 ' + '-'*40)
    test(filename='mesh_C4.nc', points=[(0., 0.1, 0.),(2., 0.1, 0.)])
    print('Test3 ' + '-'*40)
    test(filename='mesh_C4.nc', points=[(0., 1.1, 0.),(2., 1.1, 0.)])
    print('Test4 ' + '-'*40)
    test(filename='mesh_C4.nc', points=[(0., 1., 0.),(2*pi, 1., 0.)])
    print('Test5 ' + '-'*40)
    test(filename='mesh_C16.nc', points=[(0., 1.4, 0.),(2*pi, 1.4, 0.)])
    print('Test6 ' + '-'*40)
    test(filename='mesh_C4.nc', points=[(0., 1., 0.),(0., -1., 0.),(2*pi, -1., 0.)])
    print('Test7 ' + '-'*40)
    test(filename='mesh_C4.nc', points=[(3.9269908169872414, -0.6154797086703874, 0.0),(3.5342917352885173, -0.7458526607730738, 0.0)])
    print('Test8 ' + '-'*40)
    test(filename='mesh_C4.nc', points=[(3.5342917352885173, 0.7458526607730737, 0.0),(3.9269908169872414, 0.6154797086703873, 0.0)])
    print('Test9 ' + '-'*40)
    test(filename='mesh_C4.nc', points=[(0.7853981633974483, -1.04089353704597, 0.0), (0.3926990816987242, -0.7458526607730738, 0.0)])
    print('Test10' + '-'*40)
    test(filename='mesh_C4.nc', points=[(-0.0, -1.1780972450961724, 0.0), (0.7853981633974483, -1.04089353704597, 0.0)])
    print('Test11' + '-'*40)
    test(filename='mesh_C4.nc', points=[(0.0, 1.1780972450961724, 0.0), (0.7853981633974483, 1.04089353704597, 0.0)])
    print('Test12' + '-'*40)
    test(filename='mesh_C4.nc', points=[(3.9269908169872414, 1.04089353704597, 0.0), (3.141592653589793, 1.1780972450961724, 0.0)])

    # UM grid
    print('Test13' + '-'*40)
    testUm(filename='um10x4.nc', points=[(6.283185307179586, 1.1780972450961724, 0.0), (7.853981633974483, 1.5707963267948966, 0.0)], padding=4)
    print('Test14' + '-'*40)
    testUm(filename='um10x4.nc', points=[(7.853981633974483, 1.5707963267948966, 0.0), (4.71238898038469, 1.1780972450961724, 0.0)], padding=4)





