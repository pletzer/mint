from mint import Grid, ExtensiveFieldConverter
import mint
import numpy
from pathlib import Path
import vtk
from tempfile import TemporaryDirectory

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')
LAT_INDEX = 1


def test_grid():
    gr = Grid()
    filename = str(DATA_DIR / Path('test_create_grid.vtk'))
    gr.load(filename)
    ncells = gr.getNumberOfCells()
    print(f'ncells = {ncells}')
    assert(ncells == 2)
    efc = ExtensiveFieldConverter()
    efc.setGrid(gr)
    n = ncells * mint.NUM_EDGES_PER_QUAD
    vx = 1*numpy.ones(n, numpy.float64)
    vy = 2*numpy.ones(n, numpy.float64)
    dataEdge = efc.getEdgeData(vx, vy, placement=mint.CELL_BY_CELL_DATA)
    dataFace = efc.getFaceData(vx, vy, placement=mint.CELL_BY_CELL_DATA)
    for icell in range(ncells):
        for iedge in range(mint.NUM_EDGES_PER_QUAD):
            k = icell*mint.NUM_EDGES_PER_QUAD + iedge
            if iedge % 2 == 0:
                # horizontal
                assert(dataEdge[k] == 1.0)
            else:
                # vertical
                assert(dataEdge[k] == 2.0)


def test_ugrid():
    gr = Grid()
    gr.setFlags(1, 1)
    filename = str(DATA_DIR / Path('cs_4.nc'))
    gr.loadFromUgrid2DFile(f'{filename}$physics')
    nedges = gr.getNumberOfEdges()
    print(f'nedges = {nedges}')
    assert(nedges == 192)
    ncells = gr.getNumberOfCells()
    points = gr.getPoints()
    print(f'points shape: {points.shape}')

    efc = ExtensiveFieldConverter()
    efc.setGrid(gr)

    vx = 1*numpy.ones(nedges, numpy.float64) # m/s
    vy = 2*numpy.ones(nedges, numpy.float64) # m/s

    aRadius = 1.0 # units will be in radius of the planet
    for icell in range(ncells):
        for iedge in range(mint.NUM_EDGES_PER_QUAD):
            edgeId, edgeSign = gr.getEdgeId(icell, iedge)

            nodeIds = gr.getNodeIds(icell, iedge)
            p0, p1 = points[icell, iedge, :], points[icell, iedge, :]

            latmid = 0.5*(p0[LAT_INDEX] + p1[LAT_INDEX])
            cos_theta = numpy.cos(latmid * numpy.pi/180)

            print(f'icell={icell} iedge={iedge} edgeId={edgeId} nodeIds={nodeIds} p0={p0} p1={p1}')
            # convert velocity to deg/s
            vx[edgeId] *= (180/numpy.pi) / (aRadius*cos_theta)
            vy[edgeId] *= (180/numpy.pi) / aRadius

    dataEdge = efc.getEdgeData(vx, vy, placement=mint.UNIQUE_EDGE_DATA)
    dataFace = efc.getFaceData(vx, vy, placement=mint.UNIQUE_EDGE_DATA)



if __name__ == '__main__':

    test_grid()
    test_ugrid()
