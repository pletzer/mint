#!/usr/bin/env python
"""
Second-order (k=2) mimetic/Whitney edge element: circulation ("mean") +
slope DOFs per edge, built and reconstructed from scratch (mint's compiled
C++ only implements the k=1, lowest-order Whitney element -- this is NOT a
wrapper around mint.VectorInterp).

Background (see the accompanying conversation): a k=1 Whitney edge form
gives a CONSTANT tangential profile along its own edge (0th moment only --
the circulation). The k=2 (2nd order Nedelec, first kind) edge element
adds a second DOF per edge, the 1st moment ("slope"):

    t(s) = a + b*(2s - 1),   s in [0, 1] running along the edge

with a = circulation = int_0^1 t(s) ds, b = slope = int_0^1 t(s)*(2s-1) ds.
Two ways to get b from data are implemented (`mode` in buildCellDOFs):

  - 'exact':  b computed directly as the exact 1st-moment integral of the
              known field (only possible in a manufactured-field test like
              this one -- in a real pipeline you only have circulations).
  - 'split':  split the edge in two at its midpoint; b = 2*(c_far - c_near)
              recovers b EXACTLY when t(s) is purely linear, and is a
              consistent (curvature-contaminated) estimate otherwise -- see
              the module docstring of the two runner scripts for the
              derivation and how well it tracks 'exact' as h -> 0.
  - 'order1': b = 0 (plain k=1 Whitney, kept as the baseline to compare
              against).

Reference-square convention (matches the rest of this codebase, e.g.
vector_potential_reconstruction_convergence.buildExactEdgeData):

    3-->--2
    |     |
    ^     ^
    |     |
    0-->--1

Edge i0 runs corners[i0] -> corners[(i0+1) % 4], with
`sign = 1 - 2*(i0 // 2)` so that edges 0, 2 (bottom, top) both end up
oriented in the +xi direction and edges 3, 1 (left, right) both end up
oriented in the +eta direction (i.e. "both x-edges point the same way,
both y-edges point the same way").

IMPORTANT sign subtlety: `sign` flips the CIRCULATION but must NOT be
applied to the slope. Reversing which endpoint is s=0 (physically
reversing the direction of travel) negates a but leaves b UNCHANGED --
this falls out of a short derivation (reparametrize s -> u = 1-s, note
dl/du = -dl/ds, and expand): applying `sign` to b as well is exactly the
kind of silent, order-reducing bug that doesn't crash but quietly halves
the convergence rate; buildCellDOFs gets this right, and
whitney2_convergence_cartesian.py's selfTestLinearField() is a direct
regression check for it.
"""
import numpy

NUM_EDGES_PER_QUAD = 4


def gaussLegendre01(n=16):
    """n-point Gauss-Legendre nodes/weights, mapped from [-1, 1] to [0, 1]."""
    x, w = numpy.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


def edgeIntegral(pA, pB, vectorField, weightFunc=None, n=16):
    """
    Gauss-Legendre quadrature of v.dl (optionally weighted) along the
    straight chord from pA to pB.

    :param pA, pB: (2,) endpoints, in whatever 2D coordinate system
                   vectorField expects (e.g. (x, y) or (lon_deg, lat_deg))
    :param vectorField: (..., 2) points -> (..., 2) vectors, in the SAME
                   raw coordinate basis as pA/pB (no metric correction --
                   matches mint.VectorInterp's own raw-basis convention)
    :param weightFunc: None for the plain circulation (0th moment); pass
                   `lambda s: 2*s - 1` for the slope (1st moment), s in
                   [0, 1] running from pA (s=0) to pB (s=1)
    :param n: number of quadrature points (fields here are smooth
              low-order trig functions, so n=16 is already far more than
              enough -- quadrature error is negligible next to the
              reconstruction error under study; see the runner scripts)
    """
    pA = numpy.asarray(pA, dtype=numpy.float64)
    pB = numpy.asarray(pB, dtype=numpy.float64)
    s, w = gaussLegendre01(n)
    pts = pA[None, :] + s[:, None] * (pB - pA)[None, :]
    vals = vectorField(pts)
    tangent = pB - pA
    integrand = vals @ tangent
    if weightFunc is not None:
        integrand = integrand * weightFunc(s)
    return float(numpy.sum(w * integrand))


def buildCellDOFs(corners, vectorField, mode='exact', n=16):
    """
    Circulation ("a") and slope ("b") DOFs for the 4 edges of one bilinear
    quad cell -- see module docstring for the corner/edge convention and
    the sign subtlety.

    :param corners: (4, 2) array
    :param vectorField: (..., 2) -> (..., 2) callable, raw coordinate basis
    :param mode: 'order1', 'exact', or 'split' -- see module docstring
    :param n: Gauss-Legendre points per edge (or half-edge) integral
    :returns: (a, b), each shape (4,)
    """
    corners = numpy.asarray(corners, dtype=numpy.float64)
    a = numpy.zeros(4)
    b = numpy.zeros(4)
    for i0 in range(NUM_EDGES_PER_QUAD):
        i1 = (i0 + 1) % NUM_EDGES_PER_QUAD
        sign = 1 - 2 * (i0 // 2)
        pA, pB = corners[i0], corners[i1]
        a[i0] = sign * edgeIntegral(pA, pB, vectorField, None, n)
        if mode == 'order1':
            continue
        elif mode == 'exact':
            # NOTE: no `sign` here -- see module docstring. The factor of 3
            # corrects the raw moment to the actual coefficient of (2s-1)
            # in t(s) = a + b*(2s-1): (2s-1) is NOT L2-normalized on [0,1]
            # (int_0^1 (2s-1)^2 ds = 1/3, not 1), so projecting t(s) onto
            # it and dividing by that norm gives b = 3 * int t(s)(2s-1)ds.
            # Omitting this factor was the original (wrong) version and is
            # exactly what selfTestLinearField() below catches.
            b[i0] = 3.0 * edgeIntegral(pA, pB, vectorField, lambda s: 2.0 * s - 1.0, n)
        elif mode == 'split':
            mid = 0.5 * (pA + pB)
            cNear = edgeIntegral(pA, mid, vectorField, None, n)
            cFar = edgeIntegral(mid, pB, vectorField, None, n)
            b[i0] = 2.0 * (cFar - cNear)  # NOTE: no `sign` here either
        else:
            raise ValueError(f"mode must be 'order1', 'exact' or 'split', got {mode!r}")
    return a, b


def reconstructReference(a, b, xi, eta):
    """
    Reference-square reconstruction: returns (V_xi, V_eta) such that
    v.dl = V_xi*dxi + V_eta*deta at parametric point (xi, eta) in [0,1]^2.

    a[0]/b[0] (bottom) and a[2]/b[2] (top) carry the xi-direction edges;
    a[3]/b[3] (left) and a[1]/b[1] (right) the eta-direction edges.
    xi, eta may be scalars or same-shape arrays.
    """
    Vxi = ((1.0 - eta) * (a[0] + b[0] * (2.0 * xi - 1.0))
           + eta * (a[2] + b[2] * (2.0 * xi - 1.0)))
    Veta = ((1.0 - xi) * (a[3] + b[3] * (2.0 * eta - 1.0))
            + xi * (a[1] + b[1] * (2.0 * eta - 1.0)))
    return Vxi, Veta


DEFAULT_CHORD = ((0.15, 0.10), (0.75, 0.85))
"""
Default test chord in (xi, eta), deliberately NOT the (0,0)->(1,1) corner
diagonal: the plain diagonal passes exactly through the cell centroid
(0.5, 0.5), which mint's own convergence studies elsewhere in this
codebase (vector_potential_reconstruction_convergence.py) already found
to be an ISOLATED superconvergent point for pointwise reconstruction
(1st order everywhere else, 2nd order at the centroid). A path built
symmetrically around that point inherits some of that extra cancellation
and can make the k=1 ('order1') baseline look artificially close to the
k=2 ('exact'/'split') result -- this was confirmed empirically: the
diagonal gave alpha ~= 3 for ALL THREE modes, barely distinguishing them.
An off-centre, asymmetric chord like this one does not benefit from that
cancellation and cleanly separates the k=1 and k=2 convergence rates.
"""


def chordLineIntegral(a, b, xiEta0=DEFAULT_CHORD[0], xiEta1=DEFAULT_CHORD[1], n=16):
    """
    Reconstructed line integral along the straight-in-(xi, eta) chord from
    xiEta0 to xiEta1 (both in [0, 1]^2, strictly inside the cell). Purely
    in reference coordinates -- no cell geometry/Jacobian needed.
    """
    xi0, eta0 = xiEta0
    xi1, eta1 = xiEta1
    s, w = gaussLegendre01(n)
    xi = xi0 + s * (xi1 - xi0)
    eta = eta0 + s * (eta1 - eta0)
    Vxi, Veta = reconstructReference(a, b, xi, eta)
    integrand = Vxi * (xi1 - xi0) + Veta * (eta1 - eta0)
    return float(numpy.sum(w * integrand))


def bilinearMap(corners, xi, eta):
    """corners: (4, 2); xi, eta: same-shape arrays in [0, 1]. Returns (..., 2)."""
    p0, p1, p2, p3 = corners
    xi = numpy.asarray(xi, dtype=numpy.float64)[..., None]
    eta = numpy.asarray(eta, dtype=numpy.float64)[..., None]
    return ((1 - xi) * (1 - eta) * p0 + xi * (1 - eta) * p1
             + xi * eta * p2 + (1 - xi) * eta * p3)


def bilinearPartials(corners, xi, eta):
    """dX/dxi, dX/deta of the bilinear map at (xi, eta) (same-shape arrays)."""
    p0, p1, p2, p3 = corners
    xi_ = numpy.asarray(xi, dtype=numpy.float64)[..., None]
    eta_ = numpy.asarray(eta, dtype=numpy.float64)[..., None]
    dXdxi = (1 - eta_) * (p1 - p0) + eta_ * (p2 - p3)
    dXdeta = (1 - xi_) * (p3 - p0) + xi_ * (p2 - p1)
    return dXdxi, dXdeta


def exactChordLineIntegral(corners, vectorField, xiEta0=DEFAULT_CHORD[0], xiEta1=DEFAULT_CHORD[1], n=16):
    """
    TRUE line integral of vectorField along the SAME (bilinear-map) chord
    that chordLineIntegral evaluates the reconstruction along, so the two
    are directly, apples-to-apples comparable -- this measures
    reconstruction error only, not path mismatch (see module docstring of
    the runner scripts).
    """
    corners = numpy.asarray(corners, dtype=numpy.float64)
    xi0, eta0 = xiEta0
    xi1, eta1 = xiEta1
    s, w = gaussLegendre01(n)
    xi = xi0 + s * (xi1 - xi0)
    eta = eta0 + s * (eta1 - eta0)
    pts = bilinearMap(corners, xi, eta)
    dXdxi, dXdeta = bilinearPartials(corners, xi, eta)
    dXdt = dXdxi * (xi1 - xi0) + dXdeta * (eta1 - eta0)
    vals = vectorField(pts)
    integrand = numpy.sum(vals * dXdt, axis=-1)
    return float(numpy.sum(w * integrand))
