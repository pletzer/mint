#!/usr/bin/env python
"""
Rotational (non-gradient), non-divergent test velocity fields on the unit
sphere for the trajectory / contravariant-reconstruction studies
(see trajectory_contravariant_PLAN.md). Both are written as

    v(x, t) = x cross grad psi(x, t)      (tangent and divergence-free)

so the FLUX through any arc p -> q (normal n = t cross x_hat, i.e. to the
right of the direction of travel) is exactly psi(p) - psi(q). Circulations
(edge integrals v . dl, the covariant data) have a closed form only for
the solid-body field; otherwise use nodal_star.greatCircleIntegral
(Gauss-Legendre along the great-circle arc, converges spectrally).

F1 SolidBody(omega):  v = omega x x  (psi = -omega . x)
    exact edge integral along the great-circle arc p -> q (angle theta):
        omega . (p x q) * theta / sin(theta)
    exact trajectory: rotation of x0 by |omega| t about omega.

F2 NairLauritzen(k=2.4, T=5): time-dependent deformational flow
    (Nair & Lauritzen 2010, non-divergent case with a background zonal
    rotation). In the frame rotating once per T the deformation is
    multiplied by cos(pi t / T), so every tracer returns exactly to its
    starting point at t = T.

Run this file to execute the self-tests.
"""
import sys
from pathlib import Path

import numpy

sys.path.insert(0, str(Path(__file__).absolute().parent))
from nodal_star import greatCircleIntegral, arcLength  # noqa: E402


class SolidBody:

    def __init__(self, omega):
        self.omega = numpy.asarray(omega, float)

    def velocity(self, x, t=0.):
        return numpy.cross(self.omega, x)

    def psi(self, x, t=0.):
        return -numpy.sum(self.omega * x, -1)

    def edgeIntegral(self, p, q, t=0.):
        th = arcLength(p, q)
        return numpy.sum(self.omega * numpy.cross(p, q), -1) * th / numpy.sin(th)

    def basis(self):
        """Time-independent stream functions [(psi_b, gradPsi_b)]; see coeffs()."""
        om = self.omega
        return [(lambda x: -numpy.sum(om * x, -1),
                 lambda x: -numpy.broadcast_to(om, x.shape).copy())]

    def coeffs(self, t):
        return numpy.ones((1,) + numpy.shape(t))

    def trajectory(self, x0, t):
        """Rodrigues rotation of x0 by angle |omega| t about omega."""
        w = numpy.linalg.norm(self.omega)
        k = self.omega / w
        th = w * t
        return (x0 * numpy.cos(th) + numpy.cross(k, x0) * numpy.sin(th)
                + numpy.sum(k * x0, -1, keepdims=True) * k * (1. - numpy.cos(th)))


class NairLauritzen:

    def __init__(self, k=2.4, T=5.):
        self.k = k
        self.T = T

    def _parts(self, x, t):
        s = 2. * numpy.pi * t / self.T
        # cos(lat) sin(lon - s) = y cos s - x sin s
        q = x[..., 1] * numpy.cos(s) - x[..., 0] * numpy.sin(s)
        dq = numpy.stack([-numpy.sin(s) * numpy.ones_like(q),
                          numpy.cos(s) * numpy.ones_like(q),
                          numpy.zeros_like(q)], -1)
        return q, dq, self.k * numpy.cos(numpy.pi * t / self.T), 2. * numpy.pi / self.T

    def psi(self, x, t=0.):
        q, _, a, c = self._parts(x, t)
        return a * q ** 2 - c * x[..., 2]

    def gradPsi(self, x, t=0.):
        q, dq, a, c = self._parts(x, t)
        g = 2. * (a * q)[..., None] * dq
        g[..., 2] -= c
        return g

    def velocity(self, x, t=0.):
        return numpy.cross(x, self.gradPsi(x, t))

    def basis(self):
        """
        psi(x, t) = sum_b coeffs(t)[b] psi_b(x) with the time-independent
        basis y^2, x y, x^2, z, since q^2 = cos^2 s y^2 - sin 2s xy + sin^2 s x^2.
        Lets all edge data be precomputed once per basis field.
        """
        def g(k):
            return lambda x: numpy.stack(k(x), -1)
        z0 = numpy.zeros_like
        return [(lambda x: x[..., 1] ** 2, g(lambda x: (z0(x[..., 0]), 2 * x[..., 1], z0(x[..., 0])))),
                (lambda x: x[..., 0] * x[..., 1], g(lambda x: (x[..., 1], x[..., 0], z0(x[..., 0])))),
                (lambda x: x[..., 0] ** 2, g(lambda x: (2 * x[..., 0], z0(x[..., 0]), z0(x[..., 0])))),
                (lambda x: x[..., 2], g(lambda x: (z0(x[..., 0]), z0(x[..., 0]), numpy.ones_like(x[..., 0]))))]

    def coeffs(self, t):
        s = 2. * numpy.pi * t / self.T
        a = self.k * numpy.cos(numpy.pi * t / self.T)
        return numpy.stack(numpy.broadcast_arrays(a * numpy.cos(s) ** 2, -a * numpy.sin(2 * s),
                                                  a * numpy.sin(s) ** 2, -2. * numpy.pi / self.T))


def selfTest():
    rng = numpy.random.default_rng(2)

    def rpts(n):
        p = rng.normal(size=(n, 3))
        return p / numpy.linalg.norm(p, axis=-1, keepdims=True)

    p = rpts(100)
    q = p + 0.3 * rpts(100)
    q /= numpy.linalg.norm(q, axis=-1, keepdims=True)

    sb = SolidBody([0.3, -0.5, 0.8])
    # closed-form circulation vs quadrature
    quad = greatCircleIntegral(sb.velocity, p, q, nGauss=12)
    assert numpy.allclose(quad, sb.edgeIntegral(p, q), atol=1e-14, rtol=0)
    # v = x cross grad psi
    assert numpy.allclose(sb.velocity(p), numpy.cross(p, -sb.omega[None, :]))
    # trajectory: unit norm, invariant omega.x, and d/dt matches v
    x0 = rpts(50)
    t, h = 0.7, 1e-6
    xt = sb.trajectory(x0, t)
    assert numpy.allclose(numpy.linalg.norm(xt, axis=-1), 1.)
    assert numpy.allclose(xt @ sb.omega, x0 @ sb.omega)
    dxdt = (sb.trajectory(x0, t + h) - sb.trajectory(x0, t - h)) / (2 * h)
    assert numpy.allclose(dxdt, sb.velocity(xt), atol=1e-8)

    for fld, tt in ((sb, 0.), (NairLauritzen(), 1.3)):
        # flux through arc p -> q (normal to the right of travel) = psi(p) - psi(q)
        def vn(x, fld=fld, tt=tt):
            return fld.velocity(x, tt)
        th = arcLength(p, q)[..., None]
        nrm = numpy.cross(p, q)
        nrm /= numpy.linalg.norm(nrm, axis=-1, keepdims=True)
        # n = t x x_hat = -binormal of the arc, constant along it
        flux = 0.
        tk, wk = numpy.polynomial.legendre.leggauss(16)
        for tki, wki in zip(0.5 * (tk + 1), 0.5 * wk):
            x = (numpy.sin((1 - tki) * th) * p + numpy.sin(tki * th) * q) / numpy.sin(th)
            flux = flux + wki * th[..., 0] * numpy.sum(vn(x) * (-nrm), -1)
        assert numpy.allclose(flux, fld.psi(p, tt) - fld.psi(q, tt), atol=1e-13), type(fld)
        # tangent
        assert numpy.allclose(numpy.sum(vn(p) * p, -1), 0.)

    # basis decomposition reproduces psi and v
    for fld, tt in ((sb, 0.), (NairLauritzen(), 1.3)):
        c = fld.coeffs(tt)
        ps = sum(ck * b[0](p) for ck, b in zip(c, fld.basis()))
        vs = sum(ck * numpy.cross(p, b[1](p)) for ck, b in zip(c, fld.basis()))
        assert numpy.allclose(ps, fld.psi(p, tt)) and numpy.allclose(vs, fld.velocity(p, tt)), type(fld)

    # Nair-Lauritzen: gradPsi vs finite difference of psi (test only)
    nl = NairLauritzen()
    h = 1e-6
    for k in range(3):
        e = numpy.zeros(3)
        e[k] = h
        fd = (nl.psi(p + e, 1.3) - nl.psi(p - e, 1.3)) / (2 * h)
        assert numpy.allclose(fd, nl.gradPsi(p, 1.3)[:, k], atol=1e-8)
    print('trajectory_fields self-test passed')


if __name__ == '__main__':
    selfTest()
