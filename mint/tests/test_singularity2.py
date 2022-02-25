import mint
import numpy
import pytest
import vtk
from pathlib import Path

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def streamFunction(p):
    # singularity is at (0, 0)
    x, y = p[:2]
    angle = numpy.arctan2(y, x)
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

    def __init__(self, nx, ny):



        # the singularity must fall inside a cell
        assert(nx % 2 == 1)
        assert(ny % 2 == 1)

        self.nx, self.ny = nx, ny
        self.xymin = (-180., -90.)
        self.xymax = (+180., +90.)

        # create grid
        self.grid = Grid()
        # cubed-sphere
        self.grid.setFlags(1, 1)
        self.grid.loadFromUgrid2DFile(f'{DATA_DIR}/lfric_diag_wind.nc$Mesh2d')
        self.points = self.grid.getPoints()

        ncells = self.grid.getNumberOfCells()

        dx = (self.xymax[0] - self.xymin[0]) / float(nx)
        dy = (self.xymax[1] - self.xymin[1]) / float(ny)

        self.data = numpy.zeros((ncells, 4))
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

            for i0 in range(NUM_VERTS_PER_QUAD):

                i1 = (i0 + 1) % NUM_VERTS_PER_QUAD

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

        dx = (self.xymax[0] - self.xymin[0])/float(self.nx)
        dy = (self.xymax[1] - self.xymin[1])/float(self.ny)
        h = min(dx, dy)

        results = {
            'A': {'xyz': createCircle(xycenter=(0., 0.), nt=8, radius=h),
                          'flux': float('nan'),
                          'exact': 1.0,
                         },
            'B': {'xyz': createCircle(xycenter=(0., 0.), nt=32, radius=0.9*self.xymax[1]),
                           'flux': float('nan'),
                           'exact': 1.0,
                         },
            'C': {'xyz': createLoop(xybeg=(1.*self.xymax[0], +1.5*h),
                                    xyend=(1.*self.xymax[0], -1.5*h), nt=15),
                           'flux': float('nan'),
                           'exact': streamFunction((1.*self.xymax[0], -1.5*h)) - \
                                    streamFunction((1.*self.xymax[0], +1.5*h)) + 1.,
                         },
            'D': {'xyz': createLoop(xybeg=(0.91*self.xymax[0], +0.55*self.xymax[1]),
                                    xyend=(0.91*self.xymax[0], -0.55*self.xymax[1]), nt=15),
                           'flux': float('nan'),
                           'exact': streamFunction((0.91*self.xymax[0], -0.55*self.xymax[1])) - \
                                    streamFunction((0.91*self.xymax[0], +0.55*self.xymax[1])) + 1.,
                         },
        }

        for case in results:

            pli = PolylineIntegral()
            pli.setGrid(self.grid)
            # no periodicity in x
            pli.buildLocator(numCellsPerBucket=128, periodX=0, enableFolding=False)
            pli.computeWeights(results[case]['xyz'])

            results[case]['flux'] = pli.getIntegral(self.data)

            # save the contour in VTK file
            saveLineVTK(results[case]['xyz'], case + '.vtk')

            print(f'{case}: flux = {results[case]["flux"]} exact = {results[case]["exact"]} error = {results[case]["flux"] - results[case]["exact"]:.3g}')

        # for case in 'A', 'B', 'C':
        #     assert(abs(results[case]['flux'] - results[case]['exact']) < 1.e-10)
        # assert(abs(results['D']['flux'] - results['D']['exact']) < 0.03)



    def saveVectorField(self):

        nxv, nyv = 101, 101
        vi = VectorInterp()
        vi.setGrid(self.grid)
        vi.buildLocator(numCellsPerBucket=10, periodX=0.)
        xx, yy = numpy.meshgrid(numpy.linspace(self.xymin[0], self.xymax[0], nxv), 
                                numpy.linspace(self.xymin[1], self.xymax[1], nyv))
        xyz = numpy.zeros((nxv*nyv, 3), numpy.float64)
        xyz[:, 0] = xx.flat
        xyz[:, 1] = yy.flat
        vi.findPoints(xyz)
        vectors = vi.getFaceVectors(self.data, placement=0)

        ptsData = vtk.vtkDoubleArray()
        ptsData.SetNumberOfComponents(3)
        ptsData.SetNumberOfTuples(nxv * nyv)
        ptsData.SetVoidArray(xyz, nxv*nyv*3, 1)

        pts = vtk.vtkPoints()
        pts.SetData(ptsData)

        # add vector field
        vecData = vtk.vtkDoubleArray()
        vecData.SetName('vector_field')
        vecData.SetNumberOfComponents(3)
        vecData.SetNumberOfTuples(nxv * nyv)
        vecData.SetVoidArray(vectors, nxv*nyv*3, 1)

        sgrid = vtk.vtkStructuredGrid()
        sgrid.SetDimensions((nxv, nyv, 1))
        sgrid.SetPoints(pts)
        sgrid.GetPointData().AddArray(vecData)

        writer = vtk.vtkStructuredGridWriter()
        writer.SetFileName('vectors.vtk')
        writer.SetInputData(sgrid)
        writer.Update()



if __name__ == '__main__':
    cf = ContourFluxes(nx=21, ny=11)
    cf.compute()
