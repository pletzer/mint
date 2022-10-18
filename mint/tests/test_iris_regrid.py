from mint import MINTScheme, NUM_VERTS_PER_QUAD
import numpy
from pathlib import Path
import netCDF4


DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def test_sgrid():
    """
    Test regridding from data stored using SGRID conventions https://sgrid.github.io/sgrid/
    """
    with netCDF4.Dataset(DATA_DIR / Path("sgrid.nc"), 'r') as nc:
        g = nc.variables['grid']
        xname, yname = g.node_coordinates.split()
        # read the vertex coordinates
        xx = nc.variables[xname][:]
        yy = nc.variables[yname][:]
        coords = (xx, yy)
        # read the (u, v) components
        uu = nc.variables['u'][:]
        vv = nc.variables['v'][:]
        uu.units = nc.variables['u'].units
        vv.units = nc.variables['v'].units
        uv = (uu, vv)

    # target and source grids are the same
    src_coords = coords
    tgt_coords = coords

    # create a regridder
    regridder = MINTScheme(numCellsPerBucket=64, src_grid={'flags': (0, 0, 1)}).regridder(src_coords, tgt_coords)

    # regrid
    tgt_uu, tgt_vv = regridder(uv)

    # check
    n = len(uu.ravel())
    error_u = numpy.fabs(uu - tgt_uu).sum() / n
    error_v = numpy.fabs(vv - tgt_vv).sum() / n
    print(f'error u, v: {error_u} {error_v}')
    assert(error_u < 1.e-7)
    assert(error_v < 1.e-7)


if __name__ == '__main__':

    test_sgrid()

