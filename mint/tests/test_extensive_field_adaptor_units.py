"""
Unit conventions of ExtensiveFieldAdaptor.

Cell by cell data are line integrals (W1) and fluxes (W2) on the unit sphere, with
the longitude and latitude increments in radians, as expected by
VectorInterp.getVectorsOnEdges, RegridEdges.vectorApply and
PolylineIntegral.vectorGetIntegral. Unique edge data use increments in degrees, and
W2 fluxes pointing in the positive logical direction of the cell (the sign is
flipped for edges 0 and 2). On a Cartesian grid (degrees=False) the increments are
plain coordinate differences for both placements.

fromVectorField followed by toVectorField must return the original vectors, for
both placements.
"""
from pathlib import Path

import numpy
import pytest

from mint import (Grid, ExtensiveFieldAdaptor, CELL_BY_CELL_DATA, UNIQUE_EDGE_DATA,
                  FUNC_SPACE_W1, FUNC_SPACE_W2, NUM_EDGES_PER_QUAD)

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')
DEG2RAD = numpy.pi / 180.


def make_grid(lons, lats, degrees):
    nx, ny = len(lons) - 1, len(lats) - 1
    points = numpy.zeros((nx*ny, 4, 3), numpy.float64)
    k = 0
    for j in range(ny):
        for i in range(nx):
            points[k, :, :] = [(lons[i], lats[j], 0.), (lons[i + 1], lats[j], 0.),
                               (lons[i + 1], lats[j + 1], 0.), (lons[i], lats[j + 1], 0.)]
            k += 1
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=degrees)
    grid.setPoints(points)
    return grid


def from_vector_field(efa, u, v, placement, fs):
    data = numpy.zeros(u.shape)
    efa.fromVectorField(u, v, data, placement, fs)
    return data


@pytest.mark.parametrize("degrees", [True, False])
def test_cell_by_cell_round_trip(degrees):
    # away from the poles, where the W2 flux uses a different (Cartesian) estimate
    lons = numpy.linspace(0., 360., 9) if degrees else numpy.linspace(0., 2., 9)
    lats = numpy.linspace(-60., 60., 5) if degrees else numpy.linspace(-1., 1., 5)
    grid = make_grid(lons, lats, degrees)
    efa = ExtensiveFieldAdaptor()
    efa.setGrid(grid)

    n = grid.getNumberOfCells() * NUM_EDGES_PER_QUAD
    rng = numpy.random.default_rng(42)
    u, v = rng.uniform(-1., 1., n), rng.uniform(-1., 1., n)

    edge_data = from_vector_field(efa, u, v, CELL_BY_CELL_DATA, FUNC_SPACE_W1)
    face_data = from_vector_field(efa, u, v, CELL_BY_CELL_DATA, FUNC_SPACE_W2)
    u2, v2 = numpy.zeros(n), numpy.zeros(n)
    efa.toVectorField(edge_data, face_data, u2, v2, CELL_BY_CELL_DATA)

    assert numpy.abs(u2 - u).max() < 1.e-12
    assert numpy.abs(v2 - v).max() < 1.e-12


def test_cartesian_cell_by_cell_values():
    # unit square
    grid = make_grid([0., 1.], [0., 1.], degrees=False)
    efa = ExtensiveFieldAdaptor()
    efa.setGrid(grid)
    u, v = numpy.full(4, 1.0), numpy.full(4, 0.5)
    # edges go left to right and bottom to top
    edge_data = from_vector_field(efa, u, v, CELL_BY_CELL_DATA, FUNC_SPACE_W1)
    assert numpy.abs(edge_data - [1.0, 0.5, 1.0, 0.5]).max() < 1.e-14
    # flux u*dy - v*dx through each edge
    face_data = from_vector_field(efa, u, v, CELL_BY_CELL_DATA, FUNC_SPACE_W2)
    assert numpy.abs(face_data - [-0.5, 1.0, -0.5, 1.0]).max() < 1.e-14


def test_cell_by_cell_vs_unique_edge():
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=True)
    grid.loadFromUgrid2DFile(f'{DATA_DIR}/latlon4x2.nc$latlon')
    efa = ExtensiveFieldAdaptor()
    efa.setGrid(grid)
    ncells, nedges = grid.getNumberOfCells(), grid.getNumberOfEdges()
    U, V = 1.0, 0.5

    ue_w1 = from_vector_field(efa, numpy.full(nedges, U), numpy.full(nedges, V), UNIQUE_EDGE_DATA, FUNC_SPACE_W1)
    ue_w2 = from_vector_field(efa, numpy.full(nedges, U), numpy.full(nedges, V), UNIQUE_EDGE_DATA, FUNC_SPACE_W2)
    n = ncells * NUM_EDGES_PER_QUAD
    cc_w1 = from_vector_field(efa, numpy.full(n, U), numpy.full(n, V), CELL_BY_CELL_DATA, FUNC_SPACE_W1)
    cc_w2 = from_vector_field(efa, numpy.full(n, U), numpy.full(n, V), CELL_BY_CELL_DATA, FUNC_SPACE_W2)

    points = grid.getPoints()
    for icell in range(ncells):
        for ie in range(NUM_EDGES_PER_QUAD):
            k = icell*NUM_EDGES_PER_QUAD + ie
            edgeId, edgeSign = grid.getEdgeId(icell, ie)
            # radians versus degrees
            assert abs(cc_w1[k] - DEG2RAD * edgeSign * ue_w1[edgeId]) < 1.e-12
            # ... and the unique edge fluxes are flipped for edges 0 and 2, except at the
            # poles where the cell by cell flux uses a Cartesian estimate
            lat0, lat1 = points[icell, ie, 1], points[icell, (ie + 1) % 4, 1]
            if max(abs(lat0), abs(lat1)) < 90.:
                sign = -1. if ie % 2 == 0 else 1.
                assert abs(cc_w2[k] - DEG2RAD * edgeSign * sign * ue_w2[edgeId]) < 1.e-12
