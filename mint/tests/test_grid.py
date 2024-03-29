from mint import Grid
import numpy
from pathlib import Path
import vtk
from tempfile import TemporaryDirectory

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def test_create_grid():
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
                          (1., 1., 0.)]).reshape((2, 4, 3))
    gr.setPoints(points)

    # get a pointer to the points of the cell-by-cell mesh
    pts = gr.getPoints()
    print(pts)

    with TemporaryDirectory() as d:
        fname = str(Path(d) / Path('grid.vtk'))
        gr.dump(fname)


def test_attach_data():
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
                          (1., 1., 0.)]).reshape((2, 4, 3))
    gr.setPoints(points)
    # create cell data, 3 per cell
    nDataPerCell = 3
    data = numpy.arange(0, 2*nDataPerCell,
                        dtype=numpy.float64).reshape((2, nDataPerCell))
    gr.attach('mydata', data)
    with TemporaryDirectory() as d:
        fname = str(Path(d) / Path('grid.vtk'))
        gr.dump(fname)
        # read the data back to check the layout of the data
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(fname)
        reader.Update()
        ugrid = reader.GetOutput()
        arr = ugrid.GetCellData().GetArray('mydata')
        assert(arr.GetNumberOfTuples() == gr.getNumberOfCells())
        assert(arr.GetNumberOfComponents() == nDataPerCell)


def test_load_grid():
    gr = Grid()
    filename = str(DATA_DIR / Path('test_create_grid.vtk'))
    gr.load(filename)
    ncells = gr.getNumberOfCells()
    print(f'ncells = {ncells}')
    assert(ncells == 2)
    num_bad_cells = gr.check()
    assert(num_bad_cells == 0)
    points = gr.getPoints()
    assert(points.shape[0] == ncells)
    assert(points.shape[1] == 4)
    assert(points.shape[2] == 3)


def test_load_from_ugrid_file():
    gr = Grid()
    gr.setFlags(1, 1)
    filename = str(DATA_DIR / Path('cs_4.nc'))
    gr.loadFromUgrid2DFile(f'{filename}$physics')
    nedges = gr.getNumberOfEdges()
    print(f'nedges = {nedges}')
    assert(nedges == 192)
    ncells = gr.getNumberOfCells()
    for icell in range(ncells):
        for iedge in range(4):
            edgeId, edgeSign = gr.getEdgeId(icell, iedge)
            nodeIds = gr.getNodeIds(icell, iedge)
            print(f"cell {icell} edge {iedge}: " +
                  f"edgeId = {edgeId}, {edgeSign} nodeIds = {nodeIds}")
    # attaching a 3 components field to the grid
    data = numpy.array(range(ncells*4*3), numpy.float64)
    gr.attach('myData', data)
    num_bad_cells = gr.check()
    assert(num_bad_cells == 0)


def test_load_from_ugrid_file2():
    gr = Grid()
    gr.setFlags(1, 1)
    filename = str(DATA_DIR / Path('lfric_diag_wind.nc'))
    gr.loadFromUgrid2DFile(f'{filename}$Mesh2d')
    nedges = gr.getNumberOfEdges()
    print(f'nedges = {nedges}')
    assert(nedges == 3072)


def test_edge_arc_lengths():
    gr = Grid()
    gr.setFlags(1, 1)
    filename = str(DATA_DIR / Path('cs_4.nc'))
    gr.loadFromUgrid2DFile(f'{filename}$physics')
    gr.computeEdgeArcLengths()
    ncells = gr.getNumberOfCells()
    for icell in range(ncells):
        for edgeIndex in range(4):
            arcLength = gr.getEdgeArcLength(icell, edgeIndex)
            print(f""""
cell {icell} edge {edgeIndex} edge arc length/radius: {arcLength}""")


def test_load_ugrid_data():
    # a single cell
    # 3....>2....2
    # :          :
    # v          ^
    # 1          0
    # :          :
    # 0....<3....1
    xyz = numpy.array([(0.,0.,0.),
                       (1.,0.,0.),
                       (1.,1.,0.),
                       (0.,1.,0.)],
                       dtype=numpy.float64)
    face2nodes = numpy.array([(0, 1, 2, 3),], dtype=numpy.uintp)
    edge2nodes = numpy.array([(1, 2), # edge 0
                              (3, 0), # edge 1
                              (3, 2), # edge 2
                              (1, 0)],# edge 3
                              dtype=numpy.uintp)

    gr = Grid()
    gr.setFlags(0, 0)
    gr.loadFromUgrid2DData(xyz, face2nodes, edge2nodes)

    n0, n1 = gr.getNodeIds(cellId=0, edgeIndex=0)
    assert(n0 == 0)
    assert(n1 == 1)

    n0, n1 = gr.getNodeIds(cellId=0, edgeIndex=1)
    assert(n0 == 1)
    assert(n1 == 2)

    n0, n1 = gr.getNodeIds(cellId=0, edgeIndex=2)
    assert(n0 == 3)
    assert(n1 == 2)

    n0, n1 = gr.getNodeIds(cellId=0, edgeIndex=3)
    assert(n0 == 0)
    assert(n1 == 3)

    gr.dump('singleCell.vtk')


if __name__ == '__main__':

    test_edge_arc_lengths()
    test_attach_data()
    test_create_grid()
    test_load_grid()
    test_load_from_ugrid_file()
