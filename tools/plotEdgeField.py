import vtk
import numpy

vtkFileName = 'edge_upstream_00000.vtk'
varName = 'edge_integrated_velocity'
wxsize, wysize = 1800, 1000

reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(vtkFileName)
reader.Update()

grid = reader.GetOutput()
data = grid.GetCellData().GetArray(varName)
dmin, dmax = data.GetRange()
dmag = max(abs(dmin), abs(dmax))

points = grid.GetPoints()

numEdges = grid.GetNumberOfCells()

glyphs = [(vtk.vtkConeSource(), vtk.vtkPolyDataMapper(), vtk.vtkActor())  for i in range(numEdges)]

ptIds = vtk.vtkIdList()

for iEdge in range(numEdges):

    # connect
    s, m, a = glyphs[iEdge]
    m.SetInputConnection(s.GetOutputPort())
    a.SetMapper(m)

    # get the start/end points
    grid.GetCellPoints(iEdge, ptIds)
    beg = numpy.array(points.GetPoint(ptIds.GetId(0)))
    end = numpy.array(points.GetPoint(ptIds.GetId(1)))

    # radius is proportional to the intensity
    dval = data.GetTuple(iEdge)[0]
    s.SetRadius(2 * dval/dmag)

    direction = end - beg
    length = numpy.sqrt(direction.dot(direction))
    mid = 0.5*(beg + end)

    # set the cone's location
    s.SetCenter(mid)
    s.SetDirection(direction)
    s.SetHeight(length)

    # need to set color and other properties...
    # TO DO



renderer = vtk.vtkRenderer()
for g in glyphs:
    renderer.AddActor(g[2]);

renderWindow = vtk.vtkRenderWindow()
renderWindow.SetSize(wxsize, wysize)
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)


renderWindow.Render()
renderWindowInteractor.Start()

