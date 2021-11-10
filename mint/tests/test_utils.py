from mint import Grid, getIntegralsInXYZ, getIntegralsInLonLat, VectorInterp
import numpy
from pathlib import Path

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def test_pure(max_error, east, w1, xyz=True, plot=False):
    """
    Test computation of edge integrals using an eastward vector field

    :param max_error: max error vx/vy or vy/vx
    :param east: whether the vector filed is pointing east, north otherwise
    :param w1: whether the basis functions are W1 (edge) mor W2 (lateral face)
    :param xyz: whether to use Cartesian coordinates, lon-lat coordinates otherwise
    :param plot: whether to plot
    """

    gr = Grid()
    gr.setFlags(1, 1)
    filename = str(DATA_DIR / Path('lfric_diag_wind.nc'))
    gr.loadFromUgrid2D(f'{filename}$Mesh2d')
    nedge = gr.getNumberOfEdges()

    # target points
    nx, ny = 16, 8
    llon, llat = numpy.meshgrid(numpy.linspace(-170., 170., nx),
                                numpy.linspace(-80., 80., ny))
    ntarget = llon.shape[0] * llon.shape[1]
    target_points = numpy.zeros((ntarget, 3), numpy.float64)
    target_points[:, 0] = llon.flat
    target_points[:, 1] = llat.flat

    import netCDF4
    nc = netCDF4.Dataset(filename)
    # longitudes and latitudes at cell vertices
    lon = nc.variables['Mesh2d_node_x']
    lat = nc.variables['Mesh2d_node_y']
    # edge to node connectivity
    edge_node_connect = nc.variables['Mesh2d_edge_nodes']

    # only eastward velocity
    if east:
        u1 = numpy.ones((nedge,), numpy.float64)
        u2 = numpy.zeros((nedge,), numpy.float64)
    else:
        u1 = numpy.zeros((nedge,), numpy.float64)
        u2 = numpy.ones((nedge,), numpy.float64)        

    if xyz:
        # using x, y, z coords
        ue_integrated = getIntegralsInXYZ(lon, lat, edge_node_connect,
                                          u1, u2, w1=w1)
    else:
        # using multivalued lon-lat coords
        ue_integrated = getIntegralsInLonLat(lon, lat, edge_node_connect,
                                          u1, u2, w1=w1)
    # create vector interpolator
    vi = VectorInterp()
    vi.setGrid(gr)
    vi.buildLocator(numCellsPerBucket=100)
    vi.findPoints(target_points, tol2=1.e-10)

    # placement = 1 means data are on unique edges
    if w1:
        vects = vi.getEdgeVectors(ue_integrated, placement=1)
    else:
        vects = vi.getFaceVectors(ue_integrated, placement=1)

    # evaluate error
    if east:
        error = numpy.sum(numpy.fabs(vects[:, 1]/vects[:, 0]))/ntarget
    else:
        error = numpy.sum(numpy.fabs(vects[:, 0]/vects[:, 1]))/ntarget
    print(f'error = {error}')
    assert(error < max_error)

    if plot:
        from matplotlib import pylab
        pylab.quiver(target_points[:, 0], target_points[:, 1], vects[:, 0], vects[:, 1])
        pylab.title(f'east={east} w1={w1} xyz={xyz}')
        pylab.show()


if __name__ == '__main__':

    # eastward vector field
    test_pure(max_error=0.29, east=True, w1=True, xyz=True)
    test_pure(max_error=0.004, east=True, w1=False, xyz=True)
    test_pure(max_error=2.e-8, east=True, w1=True, xyz=False)
    test_pure(max_error=1.e-12, east=True, w1=False, xyz=False)

    # northward vector field
    test_pure(max_error=0.004, east=False, w1=True, xyz=True)
    test_pure(max_error=0.29, east=False, w1=False, xyz=True)
    test_pure(max_error=1.e-12, east=False, w1=True, xyz=False)
    test_pure(max_error=2.e-8, east=False, w1=False, xyz=False)