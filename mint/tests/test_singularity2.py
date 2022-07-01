import mint
import numpy
import pytest
import vtk
from pathlib import Path

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')
XCENTRE, YCENTRE = 1.9, 2.1

def streamFunction(p):
    # singularity location
    x, y = p[:2]
    angle = numpy.arctan2(y - YCENTRE, x - XCENTRE)
    return angle/(2*numpy.pi)


def createLoop(xybeg, xyend, nt):

    dt = numpy.pi/float(nt)
    a = (xyend[1] - xybeg[1])/2.
    xyz = [(xybeg[0], xybeg[1], 0.), (0., xybeg[1], 0.)] + \
          [(a*numpy.cos(t), a*numpy.sin(t), 0.) \
               for t in numpy.linspace(-numpy.pi/2. + dt, numpy.pi/2. - dt, nt-1)] + \
          [(0., xyend[1], 0.), (xyend[0], xyend[1], 0.)]

    return numpy.array(xyz)


def createCircle(xycenter, nt, radius):

    nt1 = nt + 1

    xyz = numpy.zeros((nt1, 3), numpy.float64)
    ts = numpy.linspace(0., 2*numpy.pi, nt1)
    xyz[:, 0] = xycenter[0] + radius*numpy.cos(ts)
    xyz[:, 1] = xycenter[1] + radius*numpy.sin(ts)

    return xyz


def saveLineVTK(xyz, filename):

    n = xyz.shape[0]
    ptsData = vtk.vtkDoubleArray()
    ptsData.SetNumberOfComponents(3)
    ptsData.SetNumberOfTuples(n)
    ptsData.SetVoidArray(xyz, n*3, 1)

    pts = vtk.vtkPoints()
    pts.SetData(ptsData)

    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions((n, 1, 1))
    sgrid.SetPoints(pts)

    writer = vtk.vtkStructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(sgrid)
    writer.Update()


class ContourFluxes:

    def __init__(self):


        self.xymin = (-180., -90.)
        self.xymax = (+180., +90.)

        # create grid
        self.grid = mint.Grid()
        # cubed-sphere
        self.grid.setFlags(1, 1, degrees=True)
        self.grid.loadFromUgrid2DFile(f'{DATA_DIR}/lfric_diag_wind.nc$Mesh2d')
        self.points = self.grid.getPoints()

        ncells = self.grid.getNumberOfCells()

        self.data = numpy.zeros((ncells, mint.NUM_EDGES_PER_QUAD))
        self.xymin = numpy.array([float('inf'), float('inf')])
        self.xymax = numpy.array([-float('inf'), -float('inf')])

        for k in range(ncells):

            # node indexing
            #  3-->--2
            #  |     |
            #  ^     ^
            #  |     |
            #  0-->--1

            # edge indexing
            #     2
            #  +-->--+
            #  |     |
            # 3^     ^1
            #  |     |
            #  +-->--+
            #     0

            for i0 in range(mint.NUM_VERTS_PER_QUAD):

                i1 = (i0 + 1) % mint.NUM_VERTS_PER_QUAD

                p0 = self.points[k, i0, :]
                p1 = self.points[k, i1, :]

                self.xymin[0] = min(self.xymin[0], p0[0])
                self.xymin[1] = min(self.xymin[1], p0[1])

                self.xymax[0] = max(self.xymax[0], p0[0])
                self.xymax[1] = max(self.xymax[1], p0[1])

                s0 = streamFunction(p0)
                s1 = streamFunction(p1)

                sign = 1 - 2*(i0 // 2)

                ds = sign*(s1 - s0)

                # take into account the multivaluedness
                if ds > 0.5:
                    ds -= 1.0
                elif ds < -0.5:
                    ds += 1.0

                self.data[k, i0] = ds

        self.grid.dump('grid.vtk')
        self.saveVectorField()


    def compute(self):

        # store results in dictionary
        results = {
            'A': {'xyz': createCircle(xycenter=(XCENTRE, YCENTRE), nt=8, radius=1.0),
                          'flux': float('nan'),
                          'exact': 1.0,
                         },
            'B': {'xyz': createCircle(xycenter=(0., 0.), nt=32, radius=70.0),
                           'flux': float('nan'),
                           'exact': 1.0,
                         },
            'C': {'xyz': createCircle(xycenter=(-80., -30.), nt=16, radius=30.),
                           'flux': float('nan'),
                           'exact': 0.0,
                         },
            'D': {'xyz': createLoop(xybeg=(180., +45.),
                                    xyend=(180., -45.0), nt=15),
                           'flux': float('nan'),
                           'exact': streamFunction((180., -45.)) - \
                                    streamFunction((180., +45.)) + 1.,
                         },
            'E': {'xyz': createLoop(xybeg=(162.5, +52.0),
                                    xyend=(162.5, -52.0), nt=15),
                           'flux': float('nan'),
                           'exact': streamFunction((162.5, -52.0)) - \
                                    streamFunction((162.5, +52.0)) + 1.,
                         },
        }

        targetData = []
        targetGrids = []
        resolutions = ('40x20', '80x40', '160x80')
        for res in resolutions:

            grid2 = mint.Grid()
            grid2.setFlags(0, 0, degrees=True) # uniform grid
            grid2.loadFromUgrid2DFile(f'{DATA_DIR}/latlon{res}Shifted.nc$mesh') # -180...180
            grid2.dump(f'lonlat{res}.vtk')

            regridder = mint.RegridEdges()
            regridder.setSrcGrid(self.grid)
            regridder.setDstGrid(grid2)
            regridder.buildLocator(numCellsPerBucket=100, periodX=0.0, enableFolding=0)
            regridder.computeWeights(debug=2)

            ncells2 = grid2.getNumberOfCells()
            data2 = numpy.zeros((ncells2, mint.NUM_EDGES_PER_QUAD), numpy.float64)
            regridder.apply(self.data, data2, placement=mint.CELL_BY_CELL_DATA)

            targetData.append(data2)
            targetGrids.append(grid2)

        for case in results:

            errors = []

            pli = mint.PolylineIntegral()
            pli.setGrid(self.grid)
            # no periodicity in x
            pli.buildLocator(numCellsPerBucket=128, periodX=0, enableFolding=False)
            pli.computeWeights(results[case]['xyz'])
            flux = pli.getIntegral(self.data, placement=mint.CELL_BY_CELL_DATA)

            # save the contour in VTK file
            saveLineVTK(results[case]['xyz'], case + '.vtk')

            error = flux - results[case]["exact"]
            print(f'{case} errors cs: {error:.2g}', end='')

            for ires in range(len(targetGrids)):
            
                grid2 = targetGrids[ires]
                data2 = targetData[ires]

                pli = mint.PolylineIntegral()
                pli.setGrid(grid2)
                # no periodicity in x
                pli.buildLocator(numCellsPerBucket=128, periodX=0, enableFolding=False)
                pli.computeWeights(results[case]['xyz'])
                flux2 = pli.getIntegral(data2, placement=mint.CELL_BY_CELL_DATA)

                error = flux2 - results[case]["exact"]
                print(f' {resolutions[ires]}: {error:.2g}', end='')
            print('')


    def saveVectorField(self):

        nxv, nyv = 101, 101
        vi = mint.VectorInterp()
        vi.setGrid(self.grid)
        vi.buildLocator(numCellsPerBucket=100, periodX=0.)
        xx, yy = numpy.meshgrid(numpy.linspace(self.xymin[0], self.xymax[0], nxv), 
                                numpy.linspace(self.xymin[1], self.xymax[1], nyv))
        xyz = numpy.zeros((nxv*nyv, 3), numpy.float64)
        xyz[:, 0] = xx.flat
        xyz[:, 1] = yy.flat
        vi.findPoints(xyz)
        vectors_edge = vi.getEdgeVectors(self.data, placement=mint.CELL_BY_CELL_DATA)
        vectors_face = vi.getFaceVectors(self.data, placement=mint.CELL_BY_CELL_DATA)

        ptsData = vtk.vtkDoubleArray()
        ptsData.SetNumberOfComponents(3)
        ptsData.SetNumberOfTuples(nxv * nyv)
        ptsData.SetVoidArray(xyz, nxv*nyv*3, 1)

        pts = vtk.vtkPoints()
        pts.SetData(ptsData)

        # add vector field
        vecDataEdge = vtk.vtkDoubleArray()
        vecDataEdge.SetName('vector_edge')
        vecDataEdge.SetNumberOfComponents(3)
        vecDataEdge.SetNumberOfTuples(nxv * nyv)
        vecDataEdge.SetVoidArray(vectors_edge, nxv*nyv*3, 1)

        vecDataFace = vtk.vtkDoubleArray()
        vecDataFace.SetName('vector_face')
        vecDataFace.SetNumberOfComponents(3)
        vecDataFace.SetNumberOfTuples(nxv * nyv)
        vecDataFace.SetVoidArray(vectors_face, nxv*nyv*3, 1)

        sgrid = vtk.vtkStructuredGrid()
        sgrid.SetDimensions((nxv, nyv, 1))
        sgrid.SetPoints(pts)
        sgrid.GetPointData().AddArray(vecDataEdge)
        sgrid.GetPointData().AddArray(vecDataFace)

        writer = vtk.vtkStructuredGridWriter()
        writer.SetFileName('vectors.vtk')
        writer.SetInputData(sgrid)
        writer.Update()



if __name__ == '__main__':
    cf = ContourFluxes()
    cf.compute()
    #mint.printLogMessages()

