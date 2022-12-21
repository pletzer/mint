import numpy
import iris


def saveMeshVTK(mesh, filename):
    """
    Save the mesh to a VTK file
    :param mesh: unstructurede mesh
    :param filename: file name
    """
    with open(filename, 'w') as f:
        x = mesh.node_coords.node_x.points
        y = mesh.node_coords.node_y.points
        npts = x.shape[0]
        face2node = mesh.face_node_connectivity.indices_by_location()
        ncells = face2node.shape[0]
        # write header
        f.write("# vtk DataFile Version 4.2\n")
        f.write("vtk output\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {npts} double\n")

        # write vertices
        xyz = numpy.zeros((npts, 3), numpy.float64)
        xyz[:, 0] = x
        xyz[:, 1] = y
        numpy.savetxt(f, xyz, fmt="%.10e")

        # write connectivity. Each cell is a quad made of 4 points
        cell_npts = 4
        cell_npts1 = cell_npts + 1
        f.write(f"\nCELLS {ncells} {cell_npts1*ncells}\n")
        cells = numpy.zeros((ncells, cell_npts1), numpy.int64)
        cells[:, 0] = cell_npts
        for i in range(cell_npts):
            cells[:, 1 + i] = face2node[:, i]
        numpy.savetxt(f, cells, fmt="%d")

        # write cell type
        f.write(f"\nCELL_TYPES {ncells}\n")
        # https://vtk.org/doc/release/4.2/html/vtkCellType_8h.html
        vtk_quad = 9
        cell_types = vtk_quad * numpy.ones((ncells,), numpy.int32)
        numpy.savetxt(f, cell_types, fmt="%d")



def saveVectorFieldVTK(u_cube, v_cube, filename):
    """
    Save the vectors to a VTK file
    :param u_cube: x-component
    :param v_cube: y-component
    :param filename: file name
    """
    with open(filename, 'w') as f:
        mesh = u_cube.mesh
        x = mesh.node_coords.node_x.points
        y = mesh.node_coords.node_y.points
        npts = x.shape[0]

        # write header
        f.write("# vtk DataFile Version 4.2\n")
        f.write("vtk output\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {npts} double\n")

        # write vertices
        xyz = numpy.zeros((npts, 3), numpy.float64)
        xyz[:, 0] = x
        xyz[:, 1] = y
        numpy.savetxt(f, xyz, fmt="%.10e")

        # write connectivity. Here, the cells are actually points
        cell_npts = 1
        cell_npts1 = cell_npts + 1
        f.write(f"\nCELLS {npts} {cell_npts1*npts}")
        cells = numpy.ones((npts, cell_npts1))
        cells[:, 0] = 1 # number of points defining the cell
        cells[:, 1] = range(npts)
        numpy.savetxt(f, cells, fmt="%d")

        # write cell type
        f.write(f"\nCELL_TYPES {npts}\n")
        # https://vtk.org/doc/release/4.2/html/vtkCellType_8h.html
        vtk_vertex = 1 # cell type
        cell_types = vtk_vertex * numpy.ones((npts,), numpy.int32)
        numpy.savetxt(f, cell_types, fmt="%d")



def getIntegralsInLonLat(lon, lat, edge_node_connect, u1, u2, w1=True, earth_radius=6371e3):
    """
    Convert vector components in edge/face integrals in lon-lat coordinates

    :param lon: target longitudes in degrees
    :param lat: target latitudes in degrees
    :param edge_node_connect: edge to node connectivity, array of size (nedge, 2)
    :param u1: eastward component
    :param u2: northward component
    :param w1: True if W1 (edge) field, False if W2 (face) field
    :param earth_radius: earth radius in metres
    :returns: array of size nedge
    """
    nedge = edge_node_connect.shape[0]
    # start node indices for each edge
    ibeg = edge_node_connect[:, 0]
    # end node indices for each edge
    iend = edge_node_connect[:, 1]
    # edge lengths of the edges in lon-lat
    ds = numpy.zeros((nedge, 2), numpy.float64)
    dlon = lon[iend] - lon[ibeg]
    # add/subtract a periodicity length to account for the multivaluedness
    # of the lon coordinate
    adlon = numpy.fabs(dlon)
    adlonMns = numpy.fabs(dlon - 360.)
    adlonPls = numpy.fabs(dlon + 360.)
    iMns = numpy.where(adlonMns < adlon)
    iPls = numpy.where(adlonPls < adlon)
    ds[:, 0] = dlon
    ds[iMns, 0] = dlon[iMns] - 360.
    ds[iPls, 0] = dlon[iPls] + 360.
    ds[:, 1] = lat[iend] - lat[ibeg]
    # convert ds to metres
    ds *= earth_radius * numpy.pi / 180.
    if w1:
        # integral u . dl
        return u1[:]*ds[:, 0] + u2[:]*ds[:, 1]
    else:
        # integral (zhat x u) . dl
        return u1[:]*ds[:, 1] - u2[:]*ds[:, 0]

def getIntegralsInXYZ(lon, lat, edge_node_connect, u1, u2, w1=True, earth_radius=6371e3):
    """
    Convert vector components in edge/face integrals in Cartesian coordinates

    :param lon: target longitudes in degrees
    :param lat: target latitudes in degrees
    :param edge_node_connect: edge to node connectivity, array of size (nedge, 2)
    :param u1: eastward component
    :param u2: northward component
    :param w1: True if W1 (edge) field, False if W2 (face) field
    :param earth_radius: earth radius in metres
    :returns: array of size nedge
    """
    nedge = edge_node_connect.shape[0]
    # start node indices for each edge
    ibeg = edge_node_connect[:, 0]
    # end node indices for each edge
    iend = edge_node_connect[:, 1]
    
    lonbeg = lon[ibeg]
    lonend = lon[iend]
    lonmid = (lonbeg + lonend)*0.5

    latbeg = lat[ibeg]
    latend = lat[iend]
    latmid = (latbeg + latend)*0.5

    deg2rad = numpy.pi/180.

    rho = earth_radius * numpy.cos(latbeg*deg2rad)
    xyzbeg = numpy.zeros((nedge, 3), numpy.float64)
    xyzbeg[:, 0] = rho * numpy.cos(lonbeg*deg2rad)
    xyzbeg[:, 1] = rho * numpy.sin(lonbeg*deg2rad)
    xyzbeg[:, 2] = earth_radius * numpy.sin(latbeg*deg2rad)
    
    rho = earth_radius * numpy.cos(latend*deg2rad)
    xyzend = numpy.zeros((nedge, 3), numpy.float64)
    xyzend[:, 0] = rho * numpy.cos(lonend*deg2rad)
    xyzend[:, 1] = rho * numpy.sin(lonend*deg2rad)
    xyzend[:, 2] = earth_radius * numpy.sin(latend*deg2rad)

    ds = xyzend - xyzbeg

    # unit vectors at the mid edge point
    lambda_hat = numpy.zeros((nedge, 3), numpy.float64)
    lambda_hat[:, 0] = - numpy.sin(lonmid*deg2rad)
    lambda_hat[:, 1] = numpy.cos(lonmid*deg2rad)

    theta_hat = numpy.zeros((nedge, 3), numpy.float64)
    theta_hat[:, 0] = - numpy.sin(latmid*deg2rad) * numpy.cos(lonmid*deg2rad)
    theta_hat[:, 1] = - numpy.sin(latmid*deg2rad) * numpy.sin(lonmid*deg2rad)
    theta_hat[:, 2] = numpy.cos(latmid*deg2rad)
    
    if w1:
        # integral u . dl
        return u1*numpy.sum(lambda_hat*ds, axis=1) + u2*numpy.sum(theta_hat*ds, axis=1)
    else:
        # integral (zhat x u) . dl
        return u1*numpy.sum(theta_hat*ds, axis=1) - u2*numpy.sum(lambda_hat*ds, axis=1)