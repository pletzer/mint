import numpy
import iris

DEG2RAD = numpy.pi/180.0


def computeCartesianVectors(lon, lat, u, v):
    """
    Compute the Cartesian vectors from the eastward and westward components
    :param lon: longitudes in degrees
    :param lat: latitudes in degrees 
    :param u: eastward component
    :param v: northward component
    """
    n = lon.shape[0]
    vxyz = numpy.empty((n, 3), numpy.float64)
    cosLon = numpy.cos(lon * DEG2RAD)
    sinLon = numpy.sin(lon * DEG2RAD)
    cosLat = numpy.cos(lat * DEG2RAD)
    sinLat = numpy.sin(lat * DEG2RAD)
    # vx
    vxyz[:, 0] = - u*sinLon - v*sinLat*cosLon
    # vy
    vxyz[:, 1] = + u*cosLon - v*sinLat*sinLon
    # vz
    vxyz[:, 2] = v*cosLat

    return vxyz



def computeCartesianCoords(lon, lat, radius):
    """
    Compute Cartesian coordinates from lon-lat
    :param lon: longitudes in degrees
    :param lat: latitudes in degrees
    :param radius: radius of the sphere
    :returns xyz array of size (num points, 3)
    """
    npts = lon.shape[0]
    xyz = numpy.empty((npts, 3), numpy.float64)
    cosLat = numpy.cos(lat*DEG2RAD)
    xyz[:, 0] = radius * cosLat * numpy.cos(lon*DEG2RAD)
    xyz[:, 1] = radius * cosLat * numpy.sin(lon*DEG2RAD)
    xyz[:, 2] = radius * numpy.sin(lat*DEG2RAD)
    return xyz


def computeEdgeXYZ(mesh, radius=1.0):
    """
    Compute the mid edge positions in Cartesian coordinates
    :param mesh: unstructured Iris mesh
    :param radius: radius of the sphere
    :returns set of x, y, z Cartesian coordinates
    """
    # mesh vertices in radians
    lonv_rad = mesh.node_coords.node_x.points * DEG2RAD
    latv_rad = mesh.node_coords.node_y.points * DEG2RAD

    # convert to Cartesian coords
    cosLatv = numpy.cos(latv_rad)
    xv = radius * cosLatv * numpy.cos(lonv_rad)
    yv = radius * cosLatv * numpy.sin(lonv_rad)
    zv = radius * numpy.sin(latv_rad)

    edge2node = mesh.edge_node_connectivity.indices_by_location() - \
                mesh.edge_node_connectivity.start_index

    # mid edge locations
    i0 = edge2node[:, 0]
    x0, y0, z0 = xv[i0], yv[i0], zv[i0]
    i1 = edge2node[:, 1]
    x1, y1, z1 = xv[i1], yv[i1], zv[i1]
    xe = 0.5*(x0 + x1)
    ye = 0.5*(y0 + y1)
    ze = 0.5*(z0 + z1)

    return xe, ye, ze



def computeLonLatFromXYZ(xe, ye, ze):
    """
    Compute the lon-lat coordinates from the Cartesian coordinates
    :param x: x-coordinates
    :param y: y-coordinates
    :param z: z-coordinates
    :returns set of longitudes and latitudes in degrees
    """
    rhoe = numpy.sqrt(xe*xe + ye*ye)
    lone = numpy.arctan2(ye, xe) / DEG2RAD
    late = numpy.arctan2(ze, rhoe) / DEG2RAD
    return lone, late



def saveMeshVTK(mesh, filename, radius=0.98):
    """
    Save the mesh in Cartesian coordinates to a VTK file
    :param mesh: unstructurede mesh
    :param filename: file name
    :param radius: radius of the sphere
    """
    with open(filename, 'w') as f:
        lon = mesh.node_coords.node_x.points
        lat = mesh.node_coords.node_y.points
        npts = lon.shape[0]
        # VTK uses zero based indexing
        face2node = mesh.face_node_connectivity.indices_by_location() - \
                    mesh.face_node_connectivity.start_index
        ncells = face2node.shape[0]
        # write header
        f.write("# vtk DataFile Version 4.2\n")
        f.write("vtk output\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {npts} double\n")

        # write vertices in Cartesian coordinates
        xyz = computeCartesianCoords(lon, lat, radius=radius)
        numpy.savetxt(f, xyz, fmt="%.10e")

        # write connectivity. Each cell is a quad made of 4 points
        cell_npts = 4
        cell_npts1 = cell_npts + 1
        n = cell_npts1*ncells
        f.write(f"\nCELLS {ncells} {n}\n")
        cells = numpy.empty((ncells, cell_npts1), numpy.int64)
        cells[:, 0] = cell_npts
        for i in range(cell_npts):
            cells[:, 1 + i] = face2node[:, i]
        numpy.savetxt(f, cells, fmt="%d")

        # write the cell types
        f.write(f"\nCELL_TYPES {ncells}\n")
        # https://vtk.org/doc/release/4.2/html/vtkCellType_8h.html
        vtk_quad = 9
        cell_types = vtk_quad * numpy.ones((ncells,), numpy.int32)
        numpy.savetxt(f, cell_types, fmt="%d")
        f.write("\n")



def saveVectorFieldVTK(u_cube, v_cube, filename, radius=1.0):
    """
    Save the vectors to a VTK file
    :param u_cube: eastward component
    :param v_cube: northward component
    :param filename: file name
    :param radius: radius of the sphere
    """
    with open(filename, 'w') as f:

        # extract the mesh from the cube
        mesh = u_cube.mesh

        # get the connectivity, 0-based
        edge2node = mesh.edge_node_connectivity.indices_by_location() - \
                    mesh.edge_node_connectivity.start_index

        # get the vertex coords in lon-lat coordinates
        lonv = mesh.node_coords.node_x.points
        latv = mesh.node_coords.node_y.points

        # Cartesian coordinates
        xyzv = computeCartesianCoords(lonv, latv, radius=radius)

        # compute the mid edge coords
        i0 = edge2node[:, 0]
        i1 = edge2node[:, 1]
        xyze = 0.5*(xyzv[i0, :] + xyzv[i1, :])

        # number of edges
        ne = xyze.shape[0]

        # write header
        f.write("# vtk DataFile Version 4.2\n")
        f.write("vtk output\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {ne} double\n")

        numpy.savetxt(f, xyze, fmt="%.10e")

        # write the connectivity, the cells are points
        f.write(f"\nCELLS {ne} {2*ne}\n")
        cells = numpy.empty((ne, 2))
        cells[:, 0] = 1 # number of points defining the cell
        cells[:, 1] = range(ne) # list the ne points
        numpy.savetxt(f, cells, fmt="%d")

        # write the cell types
        f.write(f"\nCELL_TYPES {ne}\n")
        # https://vtk.org/doc/release/4.2/html/vtkCellType_8h.html
        vtk_vertex = 1 # cell type
        cell_types = vtk_vertex * numpy.ones((ne,), numpy.int32)
        numpy.savetxt(f, cell_types, fmt="%d")

        # write the vectors
        f.write(f"\nPOINT_DATA {ne}\n")
        f.write(f"FIELD FieldData 1\n")
        f.write(f"vectors 3 {ne} double\n")
        u = u_cube.data
        v = v_cube.data
        xe, ye, ze = xyze[:, 0], xyze[:, 1], xyze[:, 2]
        lone = numpy.arctan2(ye, xe) / DEG2RAD
        rhoe = numpy.sqrt(xe*xe + ye*ye)
        late = numpy.arctan2(ze, rhoe) / DEG2RAD
        vxyz = computeCartesianVectors(lone, late, u_cube.data, v_cube.data)
        numpy.savetxt(f, vxyz, fmt="%.10e")

        f.write("\n")



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