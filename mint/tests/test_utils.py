from mint import Grid, getIntegralsInXYZ, getIntegralsInLonLat, VectorInterp
import numpy
import pytest
from pathlib import Path

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')
EPS = 1.22354235e-10

@pytest.fixture
def filename():
    return str(DATA_DIR / Path('lfric_diag_wind.nc'))


@pytest.fixture
def grid(filename):
    gr = Grid()
    gr.setFlags(1, 1)
    gr.loadFromUgrid2DFile(f'{filename}$Mesh2d')
    return gr


@pytest.fixture
def target_points():
    # target points
    nx, ny = 16, 8
    llon, llat = numpy.meshgrid(numpy.linspace(-170., 170., nx),
                                numpy.linspace(-80., 80., ny))
    ntarget = llon.shape[0] * llon.shape[1]
    target_points = numpy.zeros((ntarget, 3), numpy.float64)
    target_points[:, 0] = llon.flat
    target_points[:, 1] = llat.flat
    return target_points


@pytest.fixture
def connectivity(filename):
    import netCDF4
    nc = netCDF4.Dataset(filename)
    # longitudes and latitudes at cell vertices
    lon = nc.variables['Mesh2d_node_x']
    lat = nc.variables['Mesh2d_node_y']
    # edge to node connectivity
    edge_node_connect = nc.variables['Mesh2d_edge_nodes']
    return {'lon': lon, 'lat': lat, 'edge_node_connect': edge_node_connect}


def test_east_w1_xyz(grid, target_points, connectivity):
    """
    Test computation of edge integrals using an eastward vector field
    """

    nedge = grid.getNumberOfEdges()
    ntarget = target_points.shape[0]

    lon = connectivity['lon']
    lat = connectivity['lat']
    edge_node_connect = connectivity['edge_node_connect']

    # set the components
    u1 = numpy.ones((nedge,), numpy.float64)   # east
    u2 = numpy.zeros((nedge,), numpy.float64)
    ue_integrated = getIntegralsInXYZ(lon, lat, edge_node_connect,
                                      u1, u2, w1=True)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=100)
    vi.findPoints(target_points, tol2=1.e-10)

    vects = vi.getEdgeVectors(ue_integrated, placement=1)
    error = numpy.sum(numpy.fabs(vects[:, 1]/(vects[:, 0] + EPS)))/ntarget

    print(f'test_east_w1_xyz: error = {error}')
    assert(error < 0.29)

    if False:
        from matplotlib import pylab
        pylab.quiver(target_points[:, 0], target_points[:, 1],
                     vects[:, 0], vects[:, 1])
        pylab.title('test_east_w1_xyz')
        pylab.show()


def test_east_w1_lonlat(grid, target_points, connectivity):
    """
    Test computation of edge integrals using an eastward vector field
    """

    nedge = grid.getNumberOfEdges()
    ntarget = target_points.shape[0]

    lon = connectivity['lon']
    lat = connectivity['lat']
    edge_node_connect = connectivity['edge_node_connect']

    # set the components
    u1 = numpy.ones((nedge,), numpy.float64)   # east
    u2 = numpy.zeros((nedge,), numpy.float64)
    ue_integrated = getIntegralsInLonLat(lon, lat, edge_node_connect,
                                         u1, u2, w1=True)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=100)
    vi.findPoints(target_points, tol2=1.e-10)

    vects = vi.getEdgeVectors(ue_integrated, placement=1)
    error = numpy.sum(numpy.fabs(vects[:, 1]/(vects[:, 0] + EPS)))/ntarget

    print(f'test_east_w1_lonlat: error = {error}')
    assert(error < 1.e-6)

    if False:
        from matplotlib import pylab
        pylab.quiver(target_points[:, 0], target_points[:, 1],
                     vects[:, 0], vects[:, 1])
        pylab.title('test_east_w1_lonlat')
        pylab.show()


def test_east_w2_xyz(grid, target_points, connectivity):
    """
    Test computation of face integrals using an eastward vector field
    """

    nedge = grid.getNumberOfEdges()
    ntarget = target_points.shape[0]

    lon = connectivity['lon']
    lat = connectivity['lat']
    edge_node_connect = connectivity['edge_node_connect']

    # set the components
    u1 = numpy.ones((nedge,), numpy.float64)   # east
    u2 = numpy.zeros((nedge,), numpy.float64)
    ue_integrated = getIntegralsInXYZ(lon, lat, edge_node_connect,
                                      u1, u2, w1=False)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=100)
    vi.findPoints(target_points, tol2=1.e-10)

    vects = vi.getFaceVectors(ue_integrated, placement=1)
    error = numpy.sum(numpy.fabs(vects[:, 1]/(vects[:, 0] + EPS)))/ntarget

    print(f'test_east_w2_xyz: error = {error}')
    assert(error < 0.29)

    if False:
        from matplotlib import pylab
        pylab.quiver(target_points[:, 0], target_points[:, 1],
                     vects[:, 0], vects[:, 1])
        pylab.title('test_east_w2_xyz')
        pylab.show()


def test_east_w2_lonlat(grid, target_points, connectivity):
    """
    Test computation of face integrals using an eastward vector field
    """

    nedge = grid.getNumberOfEdges()
    ntarget = target_points.shape[0]

    lon = connectivity['lon']
    lat = connectivity['lat']
    edge_node_connect = connectivity['edge_node_connect']

    # set the components
    u1 = numpy.ones((nedge,), numpy.float64)   # east
    u2 = numpy.zeros((nedge,), numpy.float64)
    ue_integrated = getIntegralsInLonLat(lon, lat, edge_node_connect,
                                         u1, u2, w1=False)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=100)
    vi.findPoints(target_points, tol2=1.e-10)

    vects = vi.getFaceVectors(ue_integrated, placement=1)
    error = numpy.sum(numpy.fabs(vects[:, 1]/(vects[:, 0] + EPS)))/ntarget

    print(f'test_east_w2_lonlat: error = {error}')
    assert(error < 1.e-6)

    if False:
        from matplotlib import pylab
        pylab.quiver(target_points[:, 0], target_points[:, 1],
                     vects[:, 0], vects[:, 1])
        pylab.title('test_east_w2_lonlat')
        pylab.show()


def test_north_w1_xyz(grid, target_points, connectivity):
    """
    Test computation of edge integrals using an eastward vector field
    """

    nedge = grid.getNumberOfEdges()
    ntarget = target_points.shape[0]

    lon = connectivity['lon']
    lat = connectivity['lat']
    edge_node_connect = connectivity['edge_node_connect']

    # set the components
    u1 = numpy.zeros((nedge,), numpy.float64)
    u2 = numpy.ones((nedge,), numpy.float64)  # north
    ue_integrated = getIntegralsInXYZ(lon, lat, edge_node_connect,
                                      u1, u2, w1=True)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=100)
    vi.findPoints(target_points, tol2=1.e-10)

    vects = vi.getEdgeVectors(ue_integrated, placement=1)
    error = numpy.sum(numpy.fabs(vects[:, 0]/(vects[:, 1] + EPS)))/ntarget

    print(f'test_north_w1_xyz: error = {error}')
    assert(error < 0.004)

    if False:
        from matplotlib import pylab
        pylab.quiver(target_points[:, 0], target_points[:, 1],
                     vects[:, 0], vects[:, 1])
        pylab.title('test_north_w1_xyz')
        pylab.show()


def test_north_w1_lonlat(grid, target_points, connectivity):
    """
    Test computation of edge integrals using an eastward vector field
    """

    nedge = grid.getNumberOfEdges()
    ntarget = target_points.shape[0]

    lon = connectivity['lon']
    lat = connectivity['lat']
    edge_node_connect = connectivity['edge_node_connect']

    # set the components
    u1 = numpy.zeros((nedge,), numpy.float64)
    u2 = numpy.ones((nedge,), numpy.float64)  # north
    ue_integrated = getIntegralsInLonLat(lon, lat, edge_node_connect,
                                         u1, u2, w1=True)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=100)
    vi.findPoints(target_points, tol2=1.e-10)

    vects = vi.getEdgeVectors(ue_integrated, placement=1)
    error = numpy.sum(numpy.fabs(vects[:, 0]/(vects[:, 1] + EPS)))/ntarget

    print(f'test_north_w1_lonlat: error = {error}')
    assert(error < 1.e-6)

    if False:
        from matplotlib import pylab
        pylab.quiver(target_points[:, 0], target_points[:, 1],
                     vects[:, 0], vects[:, 1])
        pylab.title('test_north_w1_lonlat')
        pylab.show()


def test_north_w2_xyz(grid, target_points, connectivity):
    """
    Test computation of face integrals using an eastward vector field
    """

    nedge = grid.getNumberOfEdges()
    ntarget = target_points.shape[0]

    lon = connectivity['lon']
    lat = connectivity['lat']
    edge_node_connect = connectivity['edge_node_connect']

    # set the components
    u1 = numpy.zeros((nedge,), numpy.float64)
    u2 = numpy.ones((nedge,), numpy.float64)  # north
    ue_integrated = getIntegralsInXYZ(lon, lat, edge_node_connect,
                                      u1, u2, w1=False)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=100)
    vi.findPoints(target_points, tol2=1.e-10)

    vects = vi.getFaceVectors(ue_integrated, placement=1)
    error = numpy.sum(numpy.fabs(vects[:, 0]/(vects[:, 1] + EPS)))/ntarget

    print(f'test_north_w2_xyz: error = {error}')
    assert(error < 0.29)

    if False:
        from matplotlib import pylab
        pylab.quiver(target_points[:, 0], target_points[:, 1],
                     vects[:, 0], vects[:, 1])
        pylab.title('test_north_w2_xyz')
        pylab.show()


def test_north_w2_lonlat(grid, target_points, connectivity):
    """
    Test computation of face integrals using an northward vector field
    """

    nedge = grid.getNumberOfEdges()
    ntarget = target_points.shape[0]

    lon = connectivity['lon']
    lat = connectivity['lat']
    edge_node_connect = connectivity['edge_node_connect']

    # set the components
    u1 = numpy.zeros((nedge,), numpy.float64)
    u2 = numpy.ones((nedge,), numpy.float64)  # north
    ue_integrated = getIntegralsInLonLat(lon, lat, edge_node_connect,
                                         u1, u2, w1=False)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=100)
    vi.findPoints(target_points, tol2=1.e-10)

    vects = vi.getFaceVectors(ue_integrated, placement=1)
    error = numpy.sum(numpy.fabs(vects[:, 0]/(vects[:, 1] + EPS)))/ntarget

    print(f'test_north_w2_lonlat: error = {error}')
    assert(error < 1.e-6)

    if False:
        from matplotlib import pylab
        pylab.quiver(target_points[:, 0], target_points[:, 1],
                     vects[:, 0], vects[:, 1])
        pylab.title('test_north_w2_lonlat')
        pylab.show()
