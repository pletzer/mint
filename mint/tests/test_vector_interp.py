from mint import Grid
from mint import VectorInterp
from mint import CELL_BY_CELL_DATA
import numpy
import vtk
from tempfile import TemporaryDirectory
from os import sep


# noinspection SpellCheckingInspection
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


# noinspection SpellCheckingInspection
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


# noinspection SpellCheckingInspection
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


# noinspection SpellCheckingInspection
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
    vi.buildLocator(numCellsPerBucket=1, periodX=0., enableFolding=False)

    # generate targets point for the above grid
    targetPoints = generateStructuredGridPoints(nxTarget, nyTarget,
                                                v0, v1,
                                                v2, v3).reshape((-1, 3))
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

            # get the edge interpolated vectors
            vectorData = vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)
            assert(abs(vectorData.max() - 1.) < 1.e-12)
            assert(abs(vectorData.min() - 0.) < 1.e-12)

            # get the lateral flux interpolated vectors
            vectorData = vi.getFaceVectors(data, placement=CELL_BY_CELL_DATA)
            # face vectors take the sign of the area vector,
            # negative if pointing down
            assert(abs(numpy.fabs(vectorData).max() - 1.) < 1.e-12)
            assert(abs(numpy.fabs(vectorData).min() - 0.) < 1.e-12)

            # reset this edge's value back to its original
            data[cellId, edgeIndex] = 0.0


# noinspection SpellCheckingInspection
def test_rectilinear2():

    nx, ny, nxTarget, nyTarget = 1, 2, 1, 2

    v0 = numpy.array((0., 0., 0.))
    v1 = numpy.array((nx, 0., 0.))
    v2 = numpy.array((nx, ny, 0.))
    v3 = numpy.array((0., ny, 0.))
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
    vi.buildLocator(numCellsPerBucket=10, periodX=0., enableFolding=False)

    # generate targets point for the above grid
    dx = numpy.array((0.1, 0., 0.))
    dy = numpy.array((0., 0.1, 0.))
    targetPoints = generateStructuredGridPoints(nxTarget, nyTarget,
                                                v0 + dx + dy, v1 - dx + dy,
                                                v2 - dx - dy, v3 + dx - dy).reshape((-1, 3))
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

            # get the edge interpolated vectors
            vectorData = vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)

            # get the lateral flux interpolated vectors
            vectorData = vi.getFaceVectors(data, placement=CELL_BY_CELL_DATA)

            # reset this edge's value back to its original
            data[cellId, edgeIndex] = 0.0

# noinspection SpellCheckingInspection
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
    vi.buildLocator(numCellsPerBucket=1, periodX=0., enableFolding=False)

    # generate targets point for the above grid
    targetPoints = generateStructuredGridPoints(nxTarget, nyTarget,
                                                v0, v1, v2, v3).reshape((-1, 3))
    numBad = vi.findPoints(targetPoints, tol2=1.e-10)
    # all points fall within the source grid so numBad == 0
    assert(numBad == 0)

    with TemporaryDirectory() as d:
        # generate edge data
        data = numpy.zeros((numCells, 4), numpy.float64)
        for cellId in range(numCells):
            # iterate over the edges of the source grid cells
            for edgeIndex in range(4):

                # set one edge to 1, all other edges to zero
                data[cellId, edgeIndex] = 1.0

                # get the edge interpolated vectors
                vectorData = vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)
                fileName = f"{d}{sep}slanted_edgeVectors_cellId{cellId}edgeIndex{edgeIndex}Edge.vtk"
                saveVectorsVTKFile(targetPoints, vectorData, fileName)

                # get the lateral face interpolated vectors
                vectorData = vi.getFaceVectors(data, placement=CELL_BY_CELL_DATA)
                fileName = f"{d}{sep}slanted_faceVectors_cellId{cellId}edgeIndex{edgeIndex}Face.vtk"
                saveVectorsVTKFile(targetPoints, vectorData, fileName)

                # reset this edge's value back to its original
                data[cellId, edgeIndex] = 0.0


# noinspection SpellCheckingInspection
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
    vi.buildLocator(numCellsPerBucket=1, periodX=0., enableFolding=False)

    # generate targets point for the above grid
    targetPoints = generateStructuredGridPoints(nxTarget, nyTarget,
                                                v0, v1, v2, v3).reshape((-1, 3))
    numBad = vi.findPoints(targetPoints, tol2=1.e-10)
    # all points fall within the source grid so numBad == 0
    assert(numBad == 0)

    with TemporaryDirectory() as d:
        # generate edge data
        data = numpy.zeros((numCells, 4), numpy.float64)
        for cellId in range(numCells):
            # iterate over the edges of the source grid cells
            for edgeIndex in range(4):

                # set one edge to 1, all other edges to zero
                data[cellId, edgeIndex] = 1.0

                # get the edge interpolated vectors
                vectorData = vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)
                fileName = f"{d}{sep}degenerate_cellId{cellId}edgeIndex{edgeIndex}Edge.vtk"
                saveVectorsVTKFile(targetPoints, vectorData, fileName)

                # get the face interpolated vectors
                vectorData = vi.getFaceVectors(data, placement=CELL_BY_CELL_DATA)
                fileName = f"{d}{sep}degenerate_cellId{cellId}edgeIndex{edgeIndex}Face.vtk"
                saveVectorsVTKFile(targetPoints, vectorData, fileName)

                # reset this edge's value back to its original
                data[cellId, edgeIndex] = 0.0


# noinspection SpellCheckingInspection
def test_accuracy():

    #
    # internal functions
    #
    def potentialFun(lam, the):
      """
      Potential
      :param lam: lon in rad
      :param the: lat in rad
      :returns potential at lon, lat
      """
      return numpy.cos(the)*numpy.cos(lam)

    def interpPosition(vertices, xi, eta):

      return vertices[0,:]*(1-xi)*(1-eta) + \
             vertices[1,:]*xi*(1-eta) + \
             vertices[2,:]*xi*eta + \
             vertices[3,:]*(1-xi)*eta

    def computeVectorW1(vertices, edge_values, xi, eta, a=1):

      dx_dxi = (vertices[1] - vertices[0])*(1-eta) + \
               (vertices[2] - vertices[3])*eta
      dx_deta = (vertices[3] - vertices[0])*(1-xi) + \
               (vertices[2] - vertices[1])*xi

      area = numpy.cross(dx_dxi, dx_deta)
      normal = area / numpy.sqrt(area.dot(area))

      jac = area.dot(normal)

      #print(f'dx_dxi = {dx_dxi} dx_deta = {dx_deta} jac = {jac}')

      grad_xi = numpy.cross(dx_deta, normal) / jac
      grad_eta = numpy.cross(normal, dx_dxi) / jac

      # could be in any coordinate system provided
      # edge_values are consistent with the vertices coordinates
      vec = (edge_values[0]*(1-eta) + \
              edge_values[2]*eta)*grad_xi + \
             (edge_values[3]*(1-xi) + \
              edge_values[1]*xi)*grad_eta
      return vec

    def VxVyVz2UV(vx, vy, vz, lam, the):
      """
      Convert from Cartesian components to zonal-meridional
      :param vx: x component
      :param vy: y component
      :param vz: z component
      :returns u, v
      """
      u = vy*numpy.cos(lam) - vx*numpy.sin(lam)
      v = -vx*numpy.sin(the)*numpy.cos(lam) - vy*numpy.sin(the)*numpy.sin(lam) + vz*numpy.cos(the)
      return u, v

    def lamThe2XYZ(lam, the, A=1):
      x = A*numpy.cos(the)*numpy.cos(lam)
      y = A*numpy.cos(the)*numpy.sin(lam)
      z = A*numpy.sin(the)
      return numpy.array([x, y, z])


    grid = Grid()

    # grid resolution
    nx, ny = 50, 20

    nx1, ny1 = nx + 1, ny + 1
    dx, dy = 360/nx, 180/ny

    lams = numpy.linspace(0., 360, nx1)*numpy.pi/180
    thes = numpy.linspace(-90, 90, ny1)*numpy.pi/180

    numPoints = nx1 * ny1
    numCells = nx * ny
    numEdges = nx*ny1 + nx1*ny

    points = numpy.empty((numPoints, 3), numpy.float64)
    face2nodes = numpy.empty((numCells, 4), numpy.uint64)
    edge2nodes = numpy.empty((numEdges, 2), numpy.uint64)

    # build the point array
    for j in range(ny1):
        for i in range(nx1):
            points[i + j*nx1, :] = lams[i]*180/numpy.pi, thes[j]*180/numpy.pi, 0.0
    # build the connectivity arrays
    k = 0
    for j0 in range(ny):
        j1 = j0 + 1
        for i0 in range(nx):
            i1 = i0 + 1
            k0 = i0 + j0*nx1
            k1 = i1 + j0*nx1
            k2 = i1 + j1*nx1
            k3 = i0 + j1*nx1
            face2nodes[k, :] = k0, k1, k2, k3
            k += 1
    k = 0
    for j in range(ny1):
        for i0 in range(nx):
            i1 = i0 + 1
            k0 = i0 + j*nx1
            k1 = i1 + j*nx1
            edge2nodes[k, :] = k0, k1
            k += 1
    for j0 in range(ny):
        j1 = j0 + 1
        for i in range(nx1):
            k0 = i + j0*nx1
            k1 = i + j1*nx1
            edge2nodes[k, :] = k0, k1
            k += 1
    grid.loadFromUgrid2DData(points, face2nodes, edge2nodes)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=10, periodX=0., enableFolding=False)

    targetPoints = numpy.empty((numCells, 3), numpy.float64)
    edgeIntegrals = numpy.empty((numCells, 4), numpy.float64)

    errorLonLat = numpy.empty((ny, nx))
    errorCart = numpy.empty((ny, nx))



    # inside cell location
    xi, eta = 0.5, 0.5

    for j0 in range(ny):

      j1 = j0 + 1

      for i0 in range(nx):

        i1 = i0 + 1

        verts = [numpy.array([lams[i0], thes[j0], 0.]),
                 numpy.array([lams[i1], thes[j0], 0.]),
                 numpy.array([lams[i1], thes[j1], 0.]),
                 numpy.array([lams[i0], thes[j1], 0.])]

        edge_values = [
          potentialFun(lams[i1], thes[j0]) - potentialFun(lams[i0], thes[j0]), # edge 0
          potentialFun(lams[i1], thes[j1]) - potentialFun(lams[i1], thes[j0]), # edge 1
          potentialFun(lams[i1], thes[j1]) - potentialFun(lams[i0], thes[j1]), # edge 2
          potentialFun(lams[i0], thes[j1]) - potentialFun(lams[i0], thes[j0]), # edge 3
        ]

        uLonLat, vLonLat = computeVectorW1(verts, edge_values, xi=xi, eta=eta)[:2]

        # target
        lamT = lams[i0]*(1-xi)*(1-eta) + lams[i1]*xi*(1-eta) + lams[i1]*xi*eta + lams[i0]*(1-xi)*eta
        theT = thes[j0]*(1-xi)*(1-eta) + thes[j0]*xi*(1-eta) + thes[j1]*xi*eta + thes[j1]*(1-xi)*eta

        targetPoints[i0 + j0*nx, :] = numpy.array([lamT, theT, 0.]) * 180/numpy.pi
        edgeIntegrals[i0 + j0*nx, :] = edge_values

        # correct for the deg
        uLonLat /= numpy.cos(theT)

        # using Cartesian coords
        vertsXYZ = [lamThe2XYZ(verts[i][0], verts[i][1]) for i in range(4)]
        #print(vertsXYZ)
        vx, vy, vz = computeVectorW1(vertsXYZ, edge_values, xi, eta)
        uCart, vCart = VxVyVz2UV(vx, vy, vz, lamT, theT)


        uExact = - numpy.cos(theT) * numpy.sin(lamT) / numpy.cos(theT) # A = 1
        vExact = - numpy.sin(theT) * numpy.cos(lamT) # A = 1

        # print(f'lam, the = {lamT:.5f}, {theT:.5f} u, v = {uLonLat:.4f}, {vLonLat:.4f} exact u, v = {uExact:.4f}, {vExact:.4f}')
        # print(f'lam, the = {lamT:.5f}, {theT:.5f} u, v = {uCart:.4f}, {vCart:.4f} exact u, v = {uExact:.4f}, {vExact:.4f}\n\n')

        errorLonLat[j0, i0] = numpy.sqrt((uLonLat - uExact)**2 + (vLonLat - vExact)**2)
        errorCart[j0, i0] = numpy.sqrt((uCart - uExact)**2 + (vCart - vExact)**2)

    # MINT
    numBad = vi.findPoints(targetPoints, tol2=1.e-10)
    assert(numBad == 0)

    lamTarget = targetPoints[:, 0] * numpy.pi/180
    theTarget = targetPoints[:, 1] * numpy.pi/180
    uvExact = numpy.zeros((numCells, 3), numpy.float64)
    uvExact[:, 0] = - numpy.sin(lamTarget)
    uvExact[:, 1] = - numpy.sin(theTarget) * numpy.cos(lamTarget) # A = 1

    uvMint = vi.getEdgeVectors(edgeIntegrals, CELL_BY_CELL_DATA)
    # adjust the units
    uvMint *= 180/numpy.pi
    uvMint[:, 0] /= numpy.cos(theTarget)
    errorMint = numpy.sqrt((uvMint[:,0] - uvExact[:,0])**2 + (uvMint[:,1] - uvExact[:,1])**2)

    print(f'uvMint = {uvMint[numCells//3, :]}')
    print(f'uvExact = {uvExact[numCells//3, :]}')

    print(f'accuracy avg/max')
    print(f'                 lon-lat: {errorLonLat.mean():.5f}/{errorLonLat.max():.5f}')
    print(f'               cartesian: {errorCart.mean():.5f}/{errorCart.max():.5f}')
    print(f'                    mint: {errorMint.mean():.5f}/{errorMint.max():.5f}')
    assert(errorLonLat.max() < 0.004)
    assert(errorCart.max() < 0.0008)
    assert(errorMint.max() < 0.4)



if __name__ == '__main__':
    test_rectilinear()
    test_rectilinear2()
    test_slanted()
    test_degenerate()
    test_accuracy()
