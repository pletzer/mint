import numpy

class CubedSphere(object):

    def __init__(self, n):
        
        self.n = n
        nSq = n*n

        self.tile2LatLon = {
            (+1, 0, 0): numpy.zeros((nSq, 2), numpy.float64),
            (-1, 0, 0): numpy.zeros((nSq, 2), numpy.float64),
            (0, +1, 0): numpy.zeros((nSq, 2), numpy.float64),
            (0, -1, 0): numpy.zeros((nSq, 2), numpy.float64),
            (0, +1, 0): numpy.zeros((nSq, 2), numpy.float64),
            (0, -1, 0): numpy.zeros((nSq, 2), numpy.float64),
        }

        self.es = [numpy.array([1, 0, 0]),
                   numpy.array([0, 1, 0]),
                   numpy.array([0, 0, 1])]

        self.xsi, self.eta = numpy.meshgrid(range(n), range(n))


    def build(self):


        ones = numpy.ones([1, 1, 1])

        for nvect, latlons in self.tile2LatLon.items():

            # direction perpendicular to tile
            nv = numpy.array(nvect)

            uAndV = []
            for e in self.es:
                x = numpy.cross(nv, e)
                if numpy.linalg.norm(x) > 0:
                    uAndV.append(x)

            uv, vv = uAndV
            du = uv / self.n
            dv = vv / self.n
            pbeg = nv - 0.5*uv - 0.5*vv

            # x, y, z coords on the box side of size 1
            xyz = numpy.zeros((self.n, self.n, 3))
            norms = numpy.zeros((self.n, self.n))
            for i in range(3):
                xyz[..., i] = pbeg[i] + self.xsi*du[i] + self.eta*dv[i]
                norms += xyz[..., i]*xyz[..., i]
            norms = numpy.sqrt(norms)

            # project onto the sphere
            for i in range(3):
                xyz[..., i] /= norms

            xyz = xyz.reshape(self.n*self.n, 3)
            rho = numpy.sqrt(xyz[:, 0]*xyz[:, 0] + xyz[:, 1]*xyz[:, 1])

            latlons[:, 0] = numpy.arctan2(xyz[:, 2], rho) * (180./numpy.pi)
            latlons[:, 1] = numpy.arctan2(xyz[:, 1], xyz[:, 0]) * (180./numpy.pi)
            

###############################################################################

def test1():
    cs = CubedSphere(1)
    cs.build()


if __name__ == '__main__':
    test1()
