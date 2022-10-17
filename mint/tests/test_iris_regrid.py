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
        # read the vetex coordinates
        xx = nc.variables[xname][:]
        yy = nc.variables[yname][:]
        coords = (xx, yy)

    # target and source grids are the same
    src_coords = coords
    tgt_coords = coords

    # create a regridder
    regridder = MINTScheme().regridder(src_coords, tgt_coords) # add some parameters



if __name__ == '__main__':

    test_sgrid()

