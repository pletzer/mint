#!/usr/bin/env python
"""
Analytic chart utilities for the equiangular cubed sphere (unit sphere),
used by the trajectory / contravariant-reconstruction studies
(see trajectory_contravariant_PLAN.md).

Each panel f is the chart (alpha, beta) in [-pi/4, pi/4]^2,

    x(alpha, beta) = s / |s|,   s = n_f + tan(alpha) u_f + tan(beta) w_f

with (n_f, u_f, w_f) = nodal_star.FACE_FRAMES[f]. The covariant basis
a_i = dx/dxi^i, the metric g_ij = a_i . a_j, its inverse g^ij and sqrt(g) are
all known in closed form (no finite differencing). With X = tan(alpha),
Y = tan(beta), delta = 1 + X^2 + Y^2:

    g_ij    = 1/(delta^2 cos^2(alpha) cos^2(beta)) [[1+X^2, -XY], [-XY, 1+Y^2]]
    g^ij    = delta cos^2(alpha) cos^2(beta)       [[1+Y^2,  XY], [ XY, 1+X^2]]
    sqrt(g) = 1/(delta^(3/2) cos^2(alpha) cos^2(beta))

These are the same for all six panels (the frames are orthonormal). The
signed area element n_hat . (a_alpha x a_beta) = +-sqrt(g) depends on the
handedness of the panel's frame, see panelOrientation().

Run this file to execute the self-tests.
"""
import sys
from pathlib import Path

import numpy

sys.path.insert(0, str(Path(__file__).absolute().parent))
from nodal_star import FACE_FRAMES, facePoint, faceJacobian  # noqa: E402

_FRAMES = numpy.array(FACE_FRAMES, float)   # (6, 3, 3): [f, (n, u, w), xyz]


def metric(alpha, beta):
    """Covariant metric g_ij, shape (..., 2, 2)."""
    X, Y = numpy.tan(alpha), numpy.tan(beta)
    delta = 1. + X ** 2 + Y ** 2
    fac = 1. / (delta ** 2 * numpy.cos(alpha) ** 2 * numpy.cos(beta) ** 2)
    return fac[..., None, None] * numpy.stack(
        [numpy.stack([1. + X ** 2, -X * Y], -1),
         numpy.stack([-X * Y, 1. + Y ** 2], -1)], -2)


def invMetric(alpha, beta):
    """Contravariant metric g^ij, shape (..., 2, 2)."""
    X, Y = numpy.tan(alpha), numpy.tan(beta)
    delta = 1. + X ** 2 + Y ** 2
    fac = delta * numpy.cos(alpha) ** 2 * numpy.cos(beta) ** 2
    return fac[..., None, None] * numpy.stack(
        [numpy.stack([1. + Y ** 2, X * Y], -1),
         numpy.stack([X * Y, 1. + X ** 2], -1)], -2)


def sqrtMetric(alpha, beta):
    X, Y = numpy.tan(alpha), numpy.tan(beta)
    delta = 1. + X ** 2 + Y ** 2
    return 1. / (delta ** 1.5 * numpy.cos(alpha) ** 2 * numpy.cos(beta) ** 2)


def panelOrientation(face):
    """+1 if (a_alpha, a_beta, x_hat) is right-handed on this panel, else -1."""
    n, u, w = _FRAMES[face]
    return int(numpy.sign(numpy.dot(numpy.cross(u, w), n)))


def basis(face, alpha, beta):
    """Covariant basis (a_alpha, a_beta), each (..., 3). Closed form."""
    return faceJacobian(face, numpy.asarray(alpha, float), numpy.asarray(beta, float))


def dualBasis(face, alpha, beta):
    """Contravariant basis a^i = g^ij a_j, each (..., 3)."""
    aa, ab = basis(face, alpha, beta)
    gi = invMetric(alpha, beta)
    return (gi[..., 0, 0, None] * aa + gi[..., 0, 1, None] * ab,
            gi[..., 1, 0, None] * aa + gi[..., 1, 1, None] * ab)


def toChart(x):
    """
    Unit-sphere point(s) x (..., 3) -> (face, alpha, beta), choosing the
    panel whose normal is closest to x. Ties on panel boundaries go to the
    lowest face index; both charts are valid there.
    """
    x = numpy.asarray(x, float)
    proj = numpy.einsum('...k,fk->...f', x, _FRAMES[:, 0, :])
    face = numpy.argmax(proj, axis=-1)
    n, u, w = (_FRAMES[face, k] for k in range(3))
    xn = numpy.sum(x * n, -1)
    alpha = numpy.arctan(numpy.sum(x * u, -1) / xn)
    beta = numpy.arctan(numpy.sum(x * w, -1) / xn)
    return face, alpha, beta


def chartPoint(face, alpha, beta):
    return facePoint(face, numpy.asarray(alpha, float), numpy.asarray(beta, float))


def contravariant(face, alpha, beta, v):
    """Contravariant components (v^alpha, v^beta) of tangent vector(s) v."""
    da, db = dualBasis(face, alpha, beta)
    return numpy.sum(v * da, -1), numpy.sum(v * db, -1)


def selfTest():
    rng = numpy.random.default_rng(1)
    a = rng.uniform(-numpy.pi / 4, numpy.pi / 4, 200)
    b = rng.uniform(-numpy.pi / 4, numpy.pi / 4, 200)
    for f in range(6):
        aa, ab = basis(f, a, b)
        # analytic metric vs metric from the closed-form Jacobian
        gJ = numpy.stack([numpy.stack([numpy.sum(aa * aa, -1), numpy.sum(aa * ab, -1)], -1),
                          numpy.stack([numpy.sum(ab * aa, -1), numpy.sum(ab * ab, -1)], -1)], -2)
        assert numpy.allclose(gJ, metric(a, b), atol=1e-14), f
        # Jacobian vs centred finite difference (test only)
        h = 1e-6
        fd = (chartPoint(f, a + h, b) - chartPoint(f, a - h, b)) / (2 * h)
        assert numpy.allclose(fd, aa, atol=1e-8), f
        # g g^-1 = I, a_i . a^j = delta_ij, sqrt(g)
        assert numpy.allclose(metric(a, b) @ invMetric(a, b), numpy.eye(2), atol=1e-13)
        da, db = dualBasis(f, a, b)
        assert numpy.allclose(numpy.sum(aa * da, -1), 1) and numpy.allclose(numpy.sum(ab * db, -1), 1)
        assert numpy.allclose(numpy.sum(aa * db, -1), 0) and numpy.allclose(numpy.sum(ab * da, -1), 0)
        x = chartPoint(f, a, b)
        sg = numpy.sum(x * numpy.cross(aa, ab), -1)
        assert numpy.allclose(sg, panelOrientation(f) * sqrtMetric(a, b), atol=1e-14), f
        # chart round trip (interior points: no tie ambiguity)
        ff, a2, b2 = toChart(x)
        assert numpy.all(ff == f) and numpy.allclose(a2, a) and numpy.allclose(b2, b), f
    # grid lines meet at 120 deg at the cube corners (3 cells share the vertex)
    c = numpy.pi / 4
    g = metric(numpy.array(c), numpy.array(c))
    ang = numpy.degrees(numpy.arccos(g[0, 1] / numpy.sqrt(g[0, 0] * g[1, 1])))
    assert abs(ang - 120.) < 1e-10, ang
    print('cs_chart self-test passed; orientations:', [panelOrientation(f) for f in range(6)])


if __name__ == '__main__':
    selfTest()
