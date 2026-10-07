#!/usr/bin/env python
"""
Nodal vector recovery from edge line integrals ("star stencils") followed by
bilinear interpolation of the nodal vectors -- a 2nd-order-pointwise
alternative to lowest-order (k=1) Whitney reconstruction, prototyped on an
equiangular cubed sphere in 3D Cartesian coordinates (unit sphere).

This is a from-scratch Python prototype (like whitney2.py), NOT a wrapper
around mint.VectorInterp. Edges are the TRUE great-circle arcs of the
equiangular grid (every grid line of an equiangular panel is a great
circle), and the edge integrals are computed by Gauss-Legendre quadrature
along those arcs -- i.e. treated as exact data.

Method, per node n
------------------
The edges incident to n are grouped into "lines": great circles through n.
For each line with unit tangent t at n we estimate c0 = v(n).t to O(h^2),
parametrising the line by arc length sigma and writing
f(sigma) = v.t(sigma) = c0 + c1*sigma + O(sigma^2):

  * 'pair': the line has an edge on both sides of n (lengths a behind,
    b ahead, integrals J- and J+ both oriented along t):

        c0 = b/(a(a+b)) J-  +  a/(b(a+b)) J+

  * 'ray': only one side is available (cube corners; the arms of a
    panel-edge node that point into the panels, since grid lines kink
    across panel boundaries). Use the first two edges along the ray
    (lengths l1, l2, integrals I1, I2), Richardson-style:

        c0 = (2 l1 + l2)/(l1 (l1+l2)) I1  -  l1/(l2 (l1+l2)) I2

    With rayOrder=1 the plain c0 = I1/l1 is used instead (O(h)), to show
    what Richardson buys at the cube corners.

Then v(n) is the least-squares solution of v.t_k = c0_k over the lines
(sum t t^T) v = sum c0 t in n's tangent plane. Node classes on the cubed
sphere: interior (2 pairs), panel edge (1 pair along the cube edge + 2 rays),
cube corner (3 rays at 120 deg).

Reconstruction inside a cell: bilinear interpolation, in the cell's
equiangular (xi, eta), of the 4 nodal vectors' Cartesian components,
projected onto the tangent plane at the target point.
"""
import numpy

FACE_FRAMES = [
    # (normal, u, w): face point = normal + tan(alpha) u + tan(beta) w
    (( 1, 0, 0), ( 0, 1, 0), ( 0, 0, 1)),
    ((-1, 0, 0), ( 0, -1, 0), ( 0, 0, 1)),
    (( 0, 1, 0), (-1, 0, 0), ( 0, 0, 1)),
    (( 0, -1, 0), ( 1, 0, 0), ( 0, 0, 1)),
    (( 0, 0, 1), ( 0, 1, 0), (-1, 0, 0)),
    (( 0, 0, -1), ( 0, 1, 0), ( 1, 0, 0)),
]


def gaussLegendre01(n):
    x, w = numpy.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.), 0.5 * w


def facePoint(face, alpha, beta):
    """Unit-sphere point(s) of panel `face` at equiangular (alpha, beta)."""
    nrm, u, w = (numpy.array(a, float) for a in FACE_FRAMES[face])
    s = (nrm + numpy.tan(alpha)[..., None] * u + numpy.tan(beta)[..., None] * w)
    return s / numpy.linalg.norm(s, axis=-1, keepdims=True)


def faceJacobian(face, alpha, beta):
    """d x / d alpha, d x / d beta for x = facePoint(face, alpha, beta)."""
    nrm, u, w = (numpy.array(a, float) for a in FACE_FRAMES[face])
    s = (nrm + numpy.tan(alpha)[..., None] * u + numpy.tan(beta)[..., None] * w)
    r = numpy.linalg.norm(s, axis=-1, keepdims=True)
    x = s / r

    def d(ds):
        return (ds - numpy.sum(x * ds, axis=-1, keepdims=True) * x) / r
    da = d((1. / numpy.cos(alpha) ** 2)[..., None] * u)
    db = d((1. / numpy.cos(beta) ** 2)[..., None] * w)
    return da, db


def arcLength(p, q):
    return numpy.arctan2(numpy.linalg.norm(numpy.cross(p, q), axis=-1),
                         numpy.sum(p * q, axis=-1))


def greatCircleIntegral(vfield, p, q, nGauss=8):
    """
    Exact (to quadrature) int v.dl along the great-circle arc p -> q.
    p, q: (..., 3) unit vectors. vfield: (..., 3) -> (..., 3).
    """
    th = arcLength(p, q)[..., None]
    sth = numpy.sin(th)
    t, wts = gaussLegendre01(nGauss)
    res = 0.
    for tk, wk in zip(t, wts):
        x = (numpy.sin((1. - tk) * th) * p + numpy.sin(tk * th) * q) / sth
        dx = th * (-numpy.cos((1. - tk) * th) * p + numpy.cos(tk * th) * q) / sth
        res = res + wk * numpy.sum(vfield(x) * dx, axis=-1)
    return res


class CubedSphere:
    """Equiangular cubed sphere with M x M cells per panel and shared nodes."""

    def __init__(self, M, jitter=0.0, seed=0):
        """
        jitter > 0 randomly displaces the interior equiangular coordinates by
        up to +-jitter*(spacing), so neighbouring edges differ by an O(1)
        ratio at every M (grid lines stay great circles). The displacement is
        antisymmetric (alpha_{M-i} = -alpha_i) so that panels with opposite
        orientations still share their border nodes.
        """
        self.M = M
        a = numpy.linspace(-numpy.pi / 4, numpy.pi / 4, M + 1)
        if jitter > 0:
            rng = numpy.random.default_rng(seed)
            d = jitter * (a[1] - a[0]) * rng.uniform(-1, 1, M + 1)
            d = 0.5 * (d - d[::-1])
            d[0] = d[-1] = 0.
            a += d
        self.alphas = a
        A, B = numpy.meshgrid(a, a, indexing='ij')
        pts = numpy.stack([facePoint(f, A, B) for f in range(6)])  # (6, M+1, M+1, 3)
        key = numpy.round(pts.reshape(-1, 3) * 1e10).astype(numpy.int64)
        _, first, inv = numpy.unique(key, axis=0, return_index=True, return_inverse=True)
        self.nodes = pts.reshape(-1, 3)[first]
        self.nodeId = inv.reshape(6, M + 1, M + 1)   # panel (i, j) -> global node
        if len(self.nodes) != 6 * M * M + 2:
            raise RuntimeError(f'non-conforming grid: {len(self.nodes)} nodes')

        # unique undirected edges
        e = []
        nid = self.nodeId
        e.append(numpy.stack([nid[:, :-1, :], nid[:, 1:, :]], -1).reshape(-1, 2))
        e.append(numpy.stack([nid[:, :, :-1], nid[:, :, 1:]], -1).reshape(-1, 2))
        e = numpy.sort(numpy.concatenate(e), axis=1)
        self.edges = numpy.unique(e, axis=0)

        nn = len(self.nodes)
        nbrs = [[] for _ in range(nn)]
        for i, j in self.edges:
            nbrs[i].append(j)
            nbrs[j].append(i)
        self.nbrs = nbrs

        # node classes from the number of panels a node belongs to
        count = numpy.zeros(nn, int)
        for f in range(6):
            numpy.add.at(count, numpy.unique(nid[f]), 1)
        self.nodeClass = numpy.where(count == 3, 'corner',
                                     numpy.where(count == 2, 'panel_edge', 'interior'))

    def edgeIntegrals(self, vfield, nGauss=8):
        """Returns a function I(p, q) = int_{p->q} v.dl for node ids p, q sharing an edge."""
        P = self.nodes[self.edges[:, 0]]
        Q = self.nodes[self.edges[:, 1]]
        vals = greatCircleIntegral(vfield, P, Q, nGauss)
        table = {(int(i), int(j)): v for (i, j), v in zip(self.edges, vals)}

        def I(p, q):
            return table[(p, q)] if p < q else -table[(q, p)]
        return I


def _tangent(n, m):
    t = m - numpy.dot(m, n) * n
    return t / numpy.linalg.norm(t)


def nodalVectors(cs, I, rayOrder=2, nodes=None, equalWeights=False, tol=1e-9):
    """
    Nodal vectors (nnodes, 3) from edge integrals I(p, q), see module docstring.
    nodes: optional subset of node ids to compute (the others are left zero).
    equalWeights: use the equal-length formulas (J- + J+)/(a+b) and
        (3 I1 - I2)/(l1+l2) regardless of the actual lengths (for testing).
    """
    X = cs.nodes
    V = numpy.zeros_like(X)
    for n in (range(len(X)) if nodes is None else nodes):
        xn = X[n]
        nb = cs.nbrs[n]
        used = set()
        ts, cs0 = [], []
        for m in nb:
            if m in used:
                continue
            N = numpy.cross(xn, X[m])
            N /= numpy.linalg.norm(N)
            # look for a collinear neighbour on the other side of n
            opp = [k for k in nb if k != m and abs(numpy.dot(N, X[k])) < tol]
            t = _tangent(xn, X[m])
            if opp:
                k = opp[0]
                used.update((m, k))
                a = arcLength(X[k], xn)
                b = arcLength(xn, X[m])
                Jm, Jp = I(k, n), I(n, m)
                if equalWeights:
                    c0 = (Jm + Jp) / (a + b)
                else:
                    c0 = b / (a * (a + b)) * Jm + a / (b * (a + b)) * Jp
            else:
                used.add(m)
                l1 = arcLength(xn, X[m])
                I1 = I(n, m)
                if rayOrder == 1:
                    c0 = I1 / l1
                else:
                    nxt = [p for p in cs.nbrs[m] if p != n and abs(numpy.dot(N, X[p])) < tol
                           and numpy.dot(X[p] - X[m], t) > 0]
                    if not nxt:
                        raise RuntimeError(f'no collinear continuation for ray {n}->{m}')
                    p = nxt[0]
                    l2 = arcLength(X[m], X[p])
                    I2 = I(m, p)
                    if equalWeights:
                        c0 = (3 * I1 - I2) / (l1 + l2)
                    else:
                        c0 = (2 * l1 + l2) / (l1 * (l1 + l2)) * I1 - l1 / (l2 * (l1 + l2)) * I2
            ts.append(t)
            cs0.append(c0)
        T = numpy.array(ts)
        # least squares in the tangent plane: (T^T T) v = T^T c0, restricted to tangent plane
        e1 = _tangent(xn, T[0] + xn)
        e2 = numpy.cross(xn, e1)
        Tp = numpy.stack([T @ e1, T @ e2], axis=1)
        v2 = numpy.linalg.lstsq(Tp, numpy.array(cs0), rcond=None)[0]
        V[n] = v2[0] * e1 + v2[1] * e2
    return V


def cellTargets(cs, xi, eta):
    """Target points (6, M, M, 3) at equiangular cell coordinates (xi, eta)."""
    a = cs.alphas
    da = numpy.diff(a)
    A = a[:-1] + xi * da
    B = a[:-1] + eta * da
    AA, BB = numpy.meshgrid(A, B, indexing='ij')
    return numpy.stack([facePoint(f, AA, BB) for f in range(6)]), AA, BB


def bilinearInterp(cs, V, xi, eta):
    """Bilinear interpolation of nodal vectors at (xi, eta) in every cell, tangent-projected."""
    nid = cs.nodeId
    v00 = V[nid[:, :-1, :-1]]
    v10 = V[nid[:, 1:, :-1]]
    v11 = V[nid[:, 1:, 1:]]
    v01 = V[nid[:, :-1, 1:]]
    v = ((1 - xi) * (1 - eta) * v00 + xi * (1 - eta) * v10
         + xi * eta * v11 + (1 - xi) * eta * v01)
    x, _, _ = cellTargets(cs, xi, eta)
    return v - numpy.sum(v * x, axis=-1, keepdims=True) * x


def whitneyInterp(cs, I, xi, eta):
    """
    Lowest-order (k=1) Whitney reconstruction at (xi, eta) in every cell,
    from the same exact great-circle edge integrals (covariant Piola map).
    """
    nid = cs.nodeId
    M = cs.M
    Iv = numpy.vectorize(I)
    n00, n10 = nid[:, :-1, :-1], nid[:, 1:, :-1]
    n11, n01 = nid[:, 1:, 1:], nid[:, :-1, 1:]
    I0 = Iv(n00, n10)   # xi-edge at eta=0
    I2 = Iv(n01, n11)   # xi-edge at eta=1
    I3 = Iv(n00, n01)   # eta-edge at xi=0
    I1 = Iv(n10, n11)   # eta-edge at xi=1
    a = (1 - eta) * I0 + eta * I2      # coefficient of d xi
    b = (1 - xi) * I3 + xi * I1        # coefficient of d eta
    _, AA, BB = cellTargets(cs, xi, eta)
    da = numpy.diff(cs.alphas)
    out = numpy.empty((6, M, M, 3))
    for f in range(6):
        jx, jy = faceJacobian(f, AA, BB)
        jx = jx * da[:, None, None]         # d x / d xi  (cell i spacing)
        jy = jy * da[None, :, None]         # d x / d eta (cell j spacing)
        g11 = numpy.sum(jx * jx, -1)
        g12 = numpy.sum(jx * jy, -1)
        g22 = numpy.sum(jy * jy, -1)
        det = g11 * g22 - g12 ** 2
        c1 = (g22 * a[f] - g12 * b[f]) / det
        c2 = (-g12 * a[f] + g11 * b[f]) / det
        out[f] = c1[..., None] * jx + c2[..., None] * jy
    return out


# ---------------------------------------------------------------------------
# manufactured tangent field: v = P grad(chi) + x cross grad(psi)
# (P = tangent projector; the second term is tangent automatically)
# ---------------------------------------------------------------------------
def _gradChi(x):
    # chi = sin(2X + Y) * cos(Z)
    X, Y, Z = x[..., 0], x[..., 1], x[..., 2]
    s = numpy.sin(2 * X + Y)
    c = numpy.cos(2 * X + Y)
    cz = numpy.cos(Z)
    return numpy.stack([2 * c * cz, c * cz, -s * numpy.sin(Z)], -1)


def _gradPsi(x):
    # psi = cos(X - 2Z) * Y
    X, Y, Z = x[..., 0], x[..., 1], x[..., 2]
    s = numpy.sin(X - 2 * Z)
    return numpy.stack([-s * Y, numpy.cos(X - 2 * Z), 2 * s * Y], -1)


def testField(x):
    gc = _gradChi(x)
    gc = gc - numpy.sum(gc * x, -1, keepdims=True) * x
    return gc + numpy.cross(x, _gradPsi(x))
