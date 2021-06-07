from mint import Grid
from mint import VectorInterp
import numpy
import vtk


def generateStructuredGridPoints(nx, ny, v0, v1, v2, v3):
    """
    Generate structured grid points
    :param nx: number of x cells
    :param ny: number of y cells
    :param v0: south west corner
    :param v1: south east corner
    :param v2: north east corner
    :param v3: north west corner
    :returns array of size (nx, ny, 3)
    """
    # parametric
    nx1 = nx + 1
    ny1 = ny + 1
    x = numpy.linspace(0., 1., nx1)
    y = numpy.linspace(0., 1., ny1)
    xx1, yy1 = numpy.meshgrid(x, y, indexing='ij')
    xx0 = 1.0 - xx1
    yy0 = 1.0 - yy1
    # structured points
    spts = numpy.zeros(list(xx0.shape) + [3], numpy.float64)
    for j in range(3):
        spts[..., j] = xx0*yy0*v0[j] + \
                       xx1*yy0*v1[j] + \
                       xx1*yy1*v2[j] + \
                       xx0*yy1*v3[j]
    return spts


def getCellByCellPoints(spts):
    """
    Get the grid cell by cell
    :param spts: array of structured grid points, size (nx+1, ny+1, 3)
    :returns array of points, size(nx*ny, 4, 3)
    """
    nx1, ny1, _ = spts.shape
    nx, ny = nx1 - 1, ny1 - 1
    res = numpy.zeros((nx*ny, 4, 3), numpy.float64)
    res[:, 0, :] = spts[:-1, :-1, :].reshape((nx*ny, 3))
    res[:, 1, :] = spts[+1:, :-1, :].reshape((nx*ny, 3))
    res[:, 2, :] = spts[+1:, +1:, :].reshape((nx*ny, 3))
    res[:, 3, :] = spts[:-1, +1:, :].reshape((nx*ny, 3))
    return res


def saveVectorsVTKFile(spts, vectors, filename):
    """
    Save vector data in VTK file
    :param spts: array of points, size (numTargetPoints, 3)
    :param vectors: vector data of size (numTargetPoints, 3)
    :param filename: file name
    """
    npts = spts.shape[0]

    vpointData = vtk.vtkDoubleArray()
    vpointData.SetNumberOfComponents(3)
    vpointData.SetNumberOfTuples(npts)
    vpointData.SetVoidArray(spts, npts*3, 1)

    vpts = vtk.vtkPoints()
    vpts.SetData(vpointData)

    sgrid = vtk.vtkUnstructuredGrid()
    sgrid.SetPoints(vpts)

    # cloud of points
    sgrid.Allocate()
    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(1)
    for i in range(npts):
        ptIds.SetId(0, i)
        sgrid.InsertNextCell(vtk.VTK_VERTEX, ptIds)

    vdata = vtk.vtkDoubleArray()
    vdata.SetNumberOfComponents(3)
    vdata.SetNumberOfTuples(npts)
    vdata.SetVoidArray(vectors, npts*3, 1)
    vdata.SetName('basis_function_vectors')
    sgrid.GetPointData().SetScalars(vdata)

    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(sgrid)
    writer.Update()


def test_rectilinear():

    nx, ny, nxTarget, nyTarget = 1, 1, 2, 3

    v0 = (0., 0., 0.)
    v1 = (1., 0., 0.)
    v2 = (1., 1., 0.)
    v3 = (0., 1., 0.)
    # create the grid
    gr = Grid()
    cellPoints = getCellByCellPoints(generateStructuredGridPoints(nx, ny,
                                                                  v0, v1,
                                                                  v2, v3))
    gr.setPoints(cellPoints)
    numCells = cellPoints.shape[0]

    # create the interpolator
    vi = VectorInterp()
    vi.setGrid(gr)
    vi.buildLocator(numCellsPerBucket=1, periodX=0.)

    # generate targets point for the above grid
    targetPoints = generateStructuredGridPoints(nxTarget, nyTarget,
                                                v0, v1,
                                                v2, v3).reshape(-1, 3)
    numBad = vi.findPoints(targetPoints, tol2=1.e-10)
    # all points fall within the source grid so numBad == 0
    assert(numBad == 0)

    # generate edge data
    data = numpy.zeros((numCells, 4), numpy.float64)
    for cellId in range(numCells):
        # iterate over the edges of the source grid cells
        for edgeIndex in range(4):

            # set one edge to 1, all other edges to zero
            data[cellId, edgeIndex] = 1.0

            # get the interpolated vectors
            vectorData = vi.getVectors(numpy.array(data))
            assert(abs(vectorData.max() - 1.) < 1.e-12)
            assert(abs(vectorData.min() - 0.) < 1.e-12)

            # reset this edge's value back to its original
            data[cellId, edgeIndex] = 0.0


def test_slanted():

    nx, ny, nxTarget, nyTarget = 1, 1, 4, 5

    v0 = (0., 0., 0.)
    v1 = (1., 0., 0.)
    v2 = (1.5, 1., 0.)
    v3 = (0.5, 1., 0.)
    # create the grid
    gr = Grid()
    cellPoints = getCellByCellPoints(generateStructuredGridPoints(nx, ny,
                                                                  v0, v1,
                                                                  v2, v3))
    gr.setPoints(cellPoints)
    numCells = cellPoints.shape[0]

    # create the interpolator
    vi = VectorInterp()
    vi.setGrid(gr)
    vi.buildLocator(numCellsPerBucket=1, periodX=0.)

    # generate targets point for the above grid
    targetPoints = generateStructuredGridPoints(nxTarget, nyTarget,
                                                v0, v1, v2, v3).reshape(-1, 3)
    numBad = vi.findPoints(targetPoints, tol2=1.e-10)
    # all points fall within the source grid so numBad == 0
    assert(numBad == 0)

    # generate edge data
    data = numpy.zeros((numCells, 4), numpy.float64)
    for cellId in range(numCells):
        # iterate over the edges of the source grid cells
        for edgeIndex in range(4):

            # set one edge to 1, all other edges to zero
            data[cellId, edgeIndex] = 1.0

            # get the interpolated vectors
            vectorData = vi.getVectors(numpy.array(data))

            # reset this edge's value back to its original
            data[cellId, edgeIndex] = 0.0

            fileName = f'slanted_cellId{cellId}edgeIndex{edgeIndex}.vtk'
            saveVectorsVTKFile(targetPoints, vectorData, fileName)


def test_degenerate():

    nx, ny, nxTarget, nyTarget = 1, 1, 12, 11

    v0 = (0., 0., 0.)
    v1 = (1., 1., 0.)
    v2 = (1., 1., 0.)  # degenerate point
    v3 = (0., 1., 0.)
    # create the grid
    gr = Grid()
    cellPoints = getCellByCellPoints(generateStructuredGridPoints(nx, ny,
                                                                  v0, v1,
                                                                  v2, v3))
    gr.setPoints(cellPoints)
    numCells = cellPoints.shape[0]

    # create the interpolator
    vi = VectorInterp()
    vi.setGrid(gr)
    vi.buildLocator(numCellsPerBucket=1, periodX=0.)

    # generate targets point for the above grid
    targetPoints = generateStructuredGridPoints(nxTarget, nyTarget,
                                                v0, v1, v2, v3).reshape(-1, 3)
    numBad = vi.findPoints(targetPoints, tol2=1.e-10)
    # all points fall within the source grid so numBad == 0
    assert(numBad == 0)

    # generate edge data
    data = numpy.zeros((numCells, 4), numpy.float64)
    for cellId in range(numCells):
        # iterate over the edges of the source grid cells
        for edgeIndex in range(4):

            # set one edge to 1, all other edges to zero
            data[cellId, edgeIndex] = 1.0

            # get the interpolated vectors
            vectorData = vi.getVectors(numpy.array(data))

            # reset this edge's value back to its original
            data[cellId, edgeIndex] = 0.0

            fileName = f'degenerate_cellId{cellId}edgeIndex{edgeIndex}.vtk'
            saveVectorsVTKFile(targetPoints, vectorData, fileName)


if __name__ == '__main__':

    test_rectilinear()
    test_slanted()
    test_degenerate()
