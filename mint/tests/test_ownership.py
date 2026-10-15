"""
Regression tests for the ownership of grids and locators handed over to
RegridEdges and VectorInterp. Grids passed to RegridEdges.setSrcGrid/setDstGrid
and locators passed to VectorInterp.setLocator are borrowed: they belong to the
caller and must not be deleted by the regridder/interpolator. Getting this wrong
crashes the interpreter (double free or use after free), so each scenario runs
in its own Python process.
"""
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

import mint

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

# the subprocesses must import the same mint package as this process (an installed
# package, or the source tree with libmint built in place), from its parent directory
MINT_PARENT_DIR = Path(mint.__file__).absolute().parent.parent


PREAMBLE = f"""
import gc
import sys
sys.path.insert(0, {str(MINT_PARENT_DIR)!r})
import numpy
from mint import Grid, RegridEdges, VectorInterp, LonLatCellLocator

DATA_DIR = {str(DATA_DIR)!r}

def load(filename):
    g = Grid()
    g.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=True)
    g.loadFromUgrid2DFile(DATA_DIR + '/' + filename)
    return g
"""

SCENARIOS = {

    # the second setSrcGrid/setDstGrid used to delete the first grid, which
    # the caller still owns
    'set_grids_twice': """
g1 = load('latlon4x2.nc$latlon')
g2 = load('latlon4x2.nc$latlon')
rg = RegridEdges()
rg.setSrcGrid(g1)
rg.setSrcGrid(g2)
rg.setDstGrid(g1)
rg.setDstGrid(g2)
assert g1.getNumberOfCells() == 8
del rg
gc.collect()
del g1, g2
gc.collect()
""",

    # loading into a grid set with setSrcGrid used to make the regridder
    # claim ownership of it, and delete it in its destructor
    'load_after_set_grid': """
g = load('latlon4x2.nc$latlon')
rg = RegridEdges()
rg.setSrcGrid(g)
rg.setSrcGridFlags(0, 0)
rg.loadSrcGrid(f'{DATA_DIR}/latlon4x2.nc$latlon')
del rg
gc.collect()
assert g.getNumberOfCells() == 8
del g
gc.collect()
""",

    # a regridder whose grids are neither set nor loaded
    'unused_regridder': """
rg = RegridEdges()
del rg
gc.collect()
""",

    # setLocator after buildLocator used to leave the VectorInterp owning,
    # and deleting, the caller's locator
    'set_locator_after_build': """
g = load('latlon4x2.nc$latlon')
vi = VectorInterp()
vi.setGrid(g)
vi.buildLocator(numCellsPerBucket=4, periodX=360.)
loc = LonLatCellLocator()
loc.setDataSet(g)
loc.setNumberOfCellsPerBucket(4)
loc.buildLocator()
vi.setLocator(loc)
assert vi.findPoints(numpy.array([[45., 10., 0.]]), tol2=1.e-10) == 0
del vi
gc.collect()
del loc
gc.collect()
""",

    # buildLocator after setLocator replaces the borrowed locator with an
    # owned one, which must then be deleted (and the borrowed one not)
    'build_locator_after_set': """
g = load('latlon4x2.nc$latlon')
loc = LonLatCellLocator()
loc.setDataSet(g)
loc.setNumberOfCellsPerBucket(4)
loc.buildLocator()
vi = VectorInterp()
vi.setGrid(g)
vi.setLocator(loc)
vi.buildLocator(numCellsPerBucket=4, periodX=360.)
vi.buildLocator(numCellsPerBucket=4, periodX=360.)
assert vi.findPoints(numpy.array([[45., 10., 0.]]), tol2=1.e-10) == 0
del vi
gc.collect()
del loc
gc.collect()
""",
}


@pytest.mark.parametrize("name", list(SCENARIOS))
def test_ownership(name):
    script = PREAMBLE + textwrap.dedent(SCENARIOS[name]) + "\nprint('OK')\n"
    res = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True)
    print(res.stdout[-2000:])
    print(res.stderr[-2000:])
    assert res.returncode == 0, f'{name} exited with code {res.returncode}'
    assert res.stdout.strip().endswith('OK')
