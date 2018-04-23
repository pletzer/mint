import numpy

def getCornerCoords(xyz, coordIndex, n0, n1):
    xx = xyz[:, coordIndex].reshape(n0, n1)
    xx0 = xx[:-1, :-1]
    xx1 = xx[1:, :-1]
    xx2 = xx[1:, 1:]
    xx3 = xx[:-1, 1:]
    return xx0, xx1, xx2, xx3


def getTriangleAreas(xx0, xx1, xx2, yy0, yy1, yy2, zz0, zz1, zz2):
    dxx10, dyy10, dzz10 = xx1 - xx0, yy1 - yy0, zz1 - zz0
    dxx20, dyy20, dzz20 = xx2 - xx0, yy2 - yy0, zz2 - zz0
    areasxx = dyy10*dzz20 - dyy20*dzz10
    areasyy = dzz10*dxx20 - dzz20*dxx10
    areaszz = dxx10*dyy20 - dxx20*dyy10
    xx = (xx0 + xx1 + xx2)/3.
    yy = (yy0 + yy1 + yy2)/3.
    zz = (zz0 + zz1 + zz2)/3.
    rr = numpy.sqrt(xx**2 + yy**2 + zz**2)
    return (areasxx*xx + areasyy*yy + areaszz*zz)/rr


def getCellAreas(xyz, n0, n1):
    """
    Compute the cell areas of a structured 2D grid
    @param xyz arry of x,y,z points
    @param n0 structured grid size along first index
    @param n1 structured grid size along second index
    @return array of areas
    """
    xx0, xx1, xx2, xx3 = getCornerCoords(xyz, 0, n0, n1) # get x for all nodes
    yy0, yy1, yy2, yy3 = getCornerCoords(xyz, 1, n0, n1) # get y for all nodes
    zz0, zz1, zz2, zz3 = getCornerCoords(xyz, 2, n0, n1) # get z for all nodes
    # break the quad in two triangles and add the areas
    #  3---2
    #  | \ |
    #  0---1
    areas = getTriangleAreas(xx0, xx1, xx3, yy0, yy1, yy3, zz0, zz1, zz3)
    areas += getTriangleAreas(xx2, xx3, xx1, yy2, yy3, yy1, zz2, zz3, zz1)
    return areas
