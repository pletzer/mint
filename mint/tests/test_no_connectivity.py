"""
Regression tests for vector field methods called on grids built with
Grid.setPoints, which have no edge connectivity: the unique edge Ids of the
cells are unknown, so a field defined on unique edges cannot be indexed.
These methods used to index the input with uninitialised edge Ids and crash;
they must raise instead. Each call runs in its own Python process.
"""
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

SRC_DIR = Path(__file__).absolute().parent.parent.parent


PREAMBLE = f"""
import sys
sys.path.insert(0, {str(SRC_DIR)!r})
import numpy
from mint import Grid, PolylineIntegral, RegridEdges, VectorInterp, FUNC_SPACE_W1, UNIQUE_EDGE_DATA

def make_grid(nx, ny):
    points = numpy.zeros((nx*ny, 4, 3), numpy.float64)
    k = 0
    for j in range(ny):
        for i in range(nx):
            x0, x1 = 360.*i/nx, 360.*(i + 1)/nx
            y0, y1 = -90. + 180.*j/ny, -90. + 180.*(j + 1)/ny
            points[k, :, :] = [(x0, y0, 0.), (x1, y0, 0.), (x1, y1, 0.), (x0, y1, 0.)]
            k += 1
    g = Grid()
    g.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0)
    g.setPoints(points)
    return g

def expect_runtime_error(func, *args):
    try:
        func(*args)
    except RuntimeError:
        print('RAISED')
        return
    print('NOT RAISED')
"""

SCENARIOS = {

    'polyline_vectorGetIntegral': """
g = make_grid(8, 4)
pli = PolylineIntegral()
pli.setGrid(g)
pli.buildLocator(numCellsPerBucket=4, periodX=360.)
pli.computeWeights(numpy.array([(10., 10., 0.), (200., 30., 0.)]))
n = 4*g.getNumberOfCells()
u = numpy.ones(n)
v = numpy.ones(n)
expect_runtime_error(pli.vectorGetIntegral, u, v, FUNC_SPACE_W1)
""",

    'regrid_vectorApply': """
src = make_grid(8, 4)
dst = make_grid(6, 3)
rg = RegridEdges()
rg.setSrcGrid(src)
rg.setDstGrid(dst)
rg.buildLocator(numCellsPerBucket=4, periodX=360.)
rg.computeWeights()
src_u = numpy.ones(4*src.getNumberOfCells())
src_v = numpy.ones(4*src.getNumberOfCells())
dst_u = numpy.zeros(4*dst.getNumberOfCells())
dst_v = numpy.zeros(4*dst.getNumberOfCells())
expect_runtime_error(rg.vectorApply, src_u, src_v, dst_u, dst_v, FUNC_SPACE_W1)
""",

    'vectorinterp_getEdgeVectors_unique_edge_data': """
g = make_grid(8, 4)
vi = VectorInterp()
vi.setGrid(g)
vi.buildLocator(numCellsPerBucket=4, periodX=360.)
vi.findPoints(numpy.array([(10., 10., 0.), (200., 30., 0.)]), tol2=1.e-10)
data = numpy.ones(4*g.getNumberOfCells())
expect_runtime_error(vi.getEdgeVectors, data, UNIQUE_EDGE_DATA)
""",
}


@pytest.mark.parametrize("name", list(SCENARIOS))
def test_no_connectivity(name):
    script = PREAMBLE + textwrap.dedent(SCENARIOS[name])
    res = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True)
    print(res.stdout[-2000:])
    print(res.stderr[-2000:])
    assert res.returncode == 0, f'{name} exited with code {res.returncode}'
    lines = [line.strip() for line in res.stdout.splitlines()]
    assert 'RAISED' in lines and 'NOT RAISED' not in lines
