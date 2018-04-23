import numpy
import vtk
import igAreas

class CubedSphere:

    def __init__(self, numCellsPerTile, radius=1.0, coords='cartesian'):


        self.appendGrids = vtk.vtkAppendFilter()
        #appendGrids.MergePointsOn()
        self.appendGrids.SetOutputPointsPrecision(1) # double

        # create tile grids
        numPointsPerTile = numCellsPerTile + 1
        us = numpy.linspace(0., 1., numPointsPerTile)
        vs = numpy.linspace(0., 1., numPointsPerTile)
        uu, vv = numpy.meshgrid(us, vs)
        # box is [0, 1] x [0, 1], let's fit a sphere inside the box
        centre = numpy.array([0.5, 0.5, 0.5])

        self.xyzList = []
        self.areaList = []
        self.tileVertList = []
        self.tilePtsList = []
        self.tileGridList = []
        self.tileAreaList = []

        # iterate over the space dimensions
        for dim0 in range(3):

            # low or high side
            for pm in range(-1, 2, 2):

                # normal vector, pointing out
                normal = numpy.zeros((3,), numpy.float64)
                normal[dim0] = pm

                # indices of the dimensions on the tile
                dim1 = (dim0 + 1) % 3
                dim2 = (dim0 + 2) % 3

                # coordinates
                xyz = numpy.zeros((numPointsPerTile*numPointsPerTile, 3), 
                                  numpy.float64)
                # grid on the box's side/tile 
                xyz[:, dim0] = (pm + 1.0)/2.
                xyz[:, dim1] = uu.flat
                xyz[:, dim2] = vv.flat
                # fix the vertex ordering so the area points outwards
                if pm > 0:
                    xyz[:, dim1] *= -1.0
                    xyz[:, dim1] += 1.0

                # project the vertices onto sphere
                for i in range(3):
                    xyz[:, i] -= centre[i]
                dist = numpy.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
                for i in range(3):
                    # normalize
                    xyz[:, i] /= dist
                    # extend to the sphere's surface
                    xyz[:, i] *= radius

                # compute the cell areas
                areas = igAreas.getCellAreas(xyz, n0=numPointsPerTile, n1=numPointsPerTile)
                tileAreas = vtk.vtkDoubleArray()
                tileAreas.SetName('cell_areas')
                tileAreas.SetNumberOfComponents(1)
                tileAreas.SetNumberOfTuples(numCellsPerTile * numCellsPerTile)
                tileAreas.SetVoidArray(areas, numCellsPerTile * numCellsPerTile, 1)

                ntot = numPointsPerTile**2

                # create the VTK unstructured grid by combining the structured grids
                tileVerts = vtk.vtkDoubleArray()
                tileVerts.SetNumberOfComponents(3)
                tileVerts.SetNumberOfTuples(ntot)

                if coords == 'spherical':
                    # convert to lat-lon in radians
                    pass
                else:
                    tileVerts.SetVoidArray(xyz, 3*ntot, 1)

                tilePts = vtk.vtkPoints()
                tilePts.SetNumberOfPoints(ntot)
                tilePts.SetData(tileVerts)

                tileGrid = vtk.vtkStructuredGrid()
                tileGrid.SetDimensions(numPointsPerTile, numPointsPerTile, 1)
                tileGrid.SetPoints(tilePts)
                tileGrid.GetCellData().SetScalars(tileAreas)

                self.appendGrids.AddInputData(tileGrid)

                self.xyzList.append(xyz)
                self.areaList.append(areas)
                self.tileVertList.append(tileVerts)
                self.tilePtsList.append(tilePts)
                self.tileGridList.append(tileGrid)
                self.tileAreaList.append(tileAreas)

        self.appendGrids.Update()
        self.grid = self.appendGrids.GetOutput() # returns unstructured grid


    def getUnstructuredGrid(self):
        return self.grid


    def save(self, filename):
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.grid)
        writer.Update()


    def show(self):

        actors = []

        gridMapper = vtk.vtkDataSetMapper()
        gridMapper.SetInputData(self.grid)

        data = self.grid.GetCellData().GetScalars()
        if data:
            lut = vtk.vtkLookupTable()
            lut.SetHueRange(0., 0.666)
            dmin, dmax = data.GetRange()
            # reset the colour scale min to map to area 0
            dmin = 0.0
            lut.SetTableRange(dmin, dmax)
            lut.Build()

            cbar = vtk.vtkScalarBarActor()
            cbar.SetLookupTable(lut)
            #actors.append(cbar)

            gridMapper.SetLookupTable(lut)
            gridMapper.SetUseLookupTableScalarRange(1)

        gridActor = vtk.vtkActor()
        gridActor.SetMapper(gridMapper)
        #gridActor.GetProperty().SetColor(93./255., 173./255., 226./255.)
        actors.append(gridActor)

        light = vtk.vtkLight()
        light.SetFocalPoint(0., 0., 0)
        light.SetPosition(1., 1., 0.)
        #light.SetSpecularColor(0.8, 0.2, 0.0)
        #light.SetDiffuseColor(0., 0.2, 0.8)
        light.SetConeAngle(0.2)
        light.SetIntensity(1.0)

        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        # add the actors to the renderer, set the background and size
        for a in actors:
            ren.AddActor(a)
        #ren.AddLight(light)
        ren.SetBackground(1, 1, 1)
        renWin.SetSize(640, 640)
        iren.Initialize()
        renWin.Render()
        iren.Start()


    def getXYZFromLambdaThetaRho(self, ltr):
        """
        Convert from lon/lat to x, y, z
        @param ltr lambda, theta and radius
        @return Cartesian coordinates
        """
        rho = ltr[2] * numpy.cos(ltr[1])
        x = rho * numpy.cos(ltr[0])
        y = rho * numpy.sin(ltr[0])
        z = ltr[2] * numpy.sin(the)
        return numpy.array([x, y, z])


    def getLambdaThetaRhoFromXYZ(self, xyz):
        """
        Convert from Cartesian to lon/lat coordinates
        @param xyz point 
        @return longitude, latitude
        """
        x, y, z = xyz
        rhoSq = x**2 + y**2 
        r = math.sqrt(rhoSq + z**2)
        rho = math.sqrt(rhoSq)
        the = math.atan2(z, rho)
        lam = math.atan2(y, x)
        return numpy.array([lam, the, r])


#############################################################################
def test():
    numCells = 8
    cs = CubedSphere(numCells)
    grid = cs.getUnstructuredGrid()
    cs.save('cs.vtk')
    cs.show()

if __name__ == '__main__':
    test()

