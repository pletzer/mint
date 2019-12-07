import netCDF4
import argparse
import ugrid_reader
import vtk
import numpy
import random

parser = argparse.ArgumentParser(description='Check the overlap between edges')
parser.add_argument('-w', dest='weights_file', default='', 
                    help='Specify the netcdf file containing the netcdf file containing weights')
parser.add_argument('-s', dest='src_grid_file', default='', 
                    help='Specify the netcdf file containing the source grid and the mesh name in the format "FILENAME:MESHNAME"')
parser.add_argument('-d', dest='dst_grid_file', default='', 
                    help='Specify the netcdf file containing the destination grid and the mesh name in the format "FILENAME:MESHNAME"')
parser.add_argument('-S', dest='src_regularization', action='store_false', help='Disable pole regularization for src grid (required for latlon grid)')
parser.add_argument('-D', dest='dst_regularization', action='store_false', help='Disable pole regularization for dst grid (required for latlon grid)')
args = parser.parse_args()



# read the weights
w = netCDF4.Dataset(args.weights_file)
edge_param_coord_beg = w.variables['edge_param_coord_beg'][:]
edge_param_coord_end = w.variables['edge_param_coord_end'][:]
dst_cell_ids = w.variables['dst_cell_ids'][:]
src_cell_ids = w.variables['src_cell_ids'][:]
dst_face_edge_ids = w.variables['dst_face_edge_ids'][:]
src_face_edge_ids = w.variables['src_face_edge_ids'][:]
weights = w.variables['weights'][:]
print('done reading the weights')


# read the src grid
src = ugrid_reader.UgridReader(filename=args.src_grid_file, regularization=args.src_regularization)
srcGrid = src.getUnstructuredGrid()
srcPoints = srcGrid.GetPoints()
print('done reading the src grid')

# read the dst grid
dst = ugrid_reader.UgridReader(filename=args.dst_grid_file, regularization=args.dst_regularization)
dstGrid = dst.getUnstructuredGrid()
dstPoints = dstGrid.GetPoints()
print('done reading the dst grid')

dstPointBeg = numpy.zeros((3,), numpy.float64)
dstPointEnd = numpy.zeros((3,), numpy.float64)
srcPointBeg = numpy.zeros((3,), numpy.float64)
srcPointEnd = numpy.zeros((3,), numpy.float64)

dstPcoordsBeg = numpy.zeros((3,), numpy.float64)
dstPcoordsEnd = numpy.zeros((3,), numpy.float64)
srcPcoordsBeg = numpy.zeros((3,), numpy.float64)
srcPcoordsEnd = numpy.zeros((3,), numpy.float64)

amat = numpy.zeros((2,2), numpy.float64)
srhs = numpy.zeros((2,), numpy.float64)

weights = numpy.zeros((8,), numpy.float64)

subId = vtk.mutable(0)

eps = 1.236436e-10

diags = {
    'srcBegPts': [],
    'srcEndPts': [],
    'srcXpts': [],
    'dstBegPts': [],
    'dstEndPts': [],
    'dstXpts': [],
    }
numBadOverlaps = 0

# iterate over all the weights
for i in range(len(weights)):

    # dst/src cell Ids
    dstCellId = dst_cell_ids[i]
    srcCellId = src_cell_ids[i]

    # dst/src edge indices
    dstEdgeIdx = dst_face_edge_ids[i]
    srcEdgeIdx = src_face_edge_ids[i]

    # dst beg/end of edge pcoords
    dstPcoordsBeg[:] = edge_param_coord_beg[dstEdgeIdx]
    dstPcoordsEnd[:] = edge_param_coord_end[dstEdgeIdx]

    # src beg/end of edge pcoords
    srcPcoordsBeg[:] = edge_param_coord_beg[srcEdgeIdx]
    srcPcoordsEnd[:] = edge_param_coord_end[srcEdgeIdx]

    dstCell = dstGrid.GetCell(dstCellId)
    srcCell = srcGrid.GetCell(srcCellId)

    # get the dst beg/end positions
    dstCell.EvaluateLocation(subId, dstPcoordsBeg, dstPointBeg, weights)
    dstCell.EvaluateLocation(subId, dstPcoordsEnd, dstPointEnd, weights)

    # get the src beg/end positions
    srcCell.EvaluateLocation(subId, srcPcoordsBeg, srcPointBeg, weights)
    srcCell.EvaluateLocation(subId, srcPcoordsEnd, srcPointEnd, weights)

    # build the matrix system
    # column 0 = dst
    # column 1 = src
    srhs[:] = dstPointBeg[:2] - srcPointBeg[:2]
    amat[:, 0] = dstPointBeg[:2] - dstPointEnd[:2]
    amat[:, 1] = srcPointEnd[:2] - srcPointBeg[:2]

    # find the intersection
    lams = numpy.linalg.solve(amat, srhs)

    # check if intersection is along the edge
    if numpy.any(lams < -eps) or numpy.any(lams > 1. + eps):
        numBadOverlaps += 1
        print(f'no overlap for i = {i} dstCellId = {dstCellId} edge={dstEdgeIdx} srcCellId = {srcCellId} edge={srcEdgeIdx} (weight={weights[i]})')
        print(f'lambdas = {lams}')
        print(f'dst points = {dstPointBeg} -> {dstPointEnd}')
        print(f'src points = {srcPointBeg} -> {srcPointEnd}')
        print('-'*80)
        diags['srcBegPts'].append(srcPointBeg)
        diags['srcEndPts'].append(srcPointEnd)
        diags['srcXpts'].append(srcPointBeg + lams[1]*(srcPointEnd - srcPointBeg)) 
        diags['dstBegPts'].append(dstPointBeg)
        diags['dstEndPts'].append(dstPointEnd)    
        diags['dstXpts'].append(dstPointBeg + lams[0]*(dstPointEnd - dstPointBeg)) 

# create pipeline
srcCones = []
srcConeMappers = []
srcConeActors = []
dstCones = []
dstConeMappers = []
dstConeActors = []
random.seed(123)
for i in range(numBadOverlaps):

    color = (random.random(), random.random(), random.random())
    print(f'i = {i} color = {color}')

    # src edge
    scone = vtk.vtkConeSource()
    sconem = vtk.vtkPolyDataMapper()
    sconea = vtk.vtkActor()
    sconem.SetInputConnection(scone.GetOutputPort())
    sconea.SetMapper(sconem)

    scone.SetCenter(0.5*(diags['srcBegPts'][i] + diags['srcEndPts'][i]))
    scone.SetDirection(diags['srcEndPts'][i] - diags['srcBegPts'][i])
    scone.SetRadius(0.2)
    sconea.GetProperty().SetColor(color)

    srcCones.append(scone)
    srcConeMappers.append(sconem)
    srcConeActors.append(sconea)

    # dst edge
    dcone = vtk.vtkConeSource()
    dconem = vtk.vtkPolyDataMapper()
    dconea = vtk.vtkActor()
    dconem.SetInputConnection(dcone.GetOutputPort())
    dconea.SetMapper(dconem)

    dcone.SetCenter(0.5*(diags['dstBegPts'][i] + diags['dstEndPts'][i]))
    dcone.SetDirection(diags['dstEndPts'][i] - diags['dstBegPts'][i])
    dcone.SetRadius(0.5)
    dconea.GetProperty().SetColor(color)

    dstCones.append(dcone)
    dstConeMappers.append(dconem)
    dstConeActors.append(dconea)

renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindowInteractor = vtk.vtkRenderWindowInteractor()

renderWindow.AddRenderer(renderer)
renderWindowInteractor.SetRenderWindow(renderWindow)

for a in srcConeActors + dstConeActors:
    renderer.AddActor(a)

renderWindow.Render()
renderWindowInteractor.Start()



