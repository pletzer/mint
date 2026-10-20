#!/usr/bin/env python
"""
Cell-by-cell trajectory integrator in the parametric (index) space of an
equiangular cubed sphere -- Phase 3 of trajectory_contravariant_PLAN.md.

A particle is (face f, cell i, j, local xi, eta in [0,1]^2). Inside a cell
it obeys d xi/dt = v^xi, d eta/dt = v^eta (cell-local contravariant
components, from one of the reconstruction arms below) and is advanced
with RK4. A step that leaves the cell is shortened so that it ends exactly
on the cell boundary (exit time from the cubic Hermite interpolant of the
step, then one RK4 step of that length), and the particle is handed to the
neighbouring cell -- across panels via xyz and the neighbour's chart.

Whitney (arms A, B) fields have discontinuous normal components, so the
normal velocity can point INTO a shared edge from both sides. Such a
particle "slides": it stays on the edge and moves with its tangential
component (which IS continuous for Whitney) until one side lets it go.

Arms (cell-local contravariant velocity):
  A  Whitney k=1 from circulations, raised with the exact metric
  B  Whitney k=1 from circulations, raised with the bilinear-corner metric
  C  RT0 from fluxes / exact Jacobian
  D  nodal star vectors, bilinear, projected onto the exact dual basis
  W  Whitney k=2 (whitney2.py edge-only element: circulation + slope per
     edge, tangential profile linear along each edge, no interior bubble
     DOFs), slopes from exact first moments in the cell coordinate, raised
     with the exact metric
  R  mint's order 2 element (mint.WhitneyElement(2, ...), Nedelec degree 2 on
     2x2 blocks of cells with Rapetti-Bossavit small-edge DOFs): the same
     circulations as A, no extra data, raised with the exact metric. Evaluated
     in the panel chart, where the equiangular grid is uniform, so the block's
     piecewise cell coordinates (biquadratic=False) are its exact geometry
  S  as R, raised with the biquadratic-block metric: the covariant basis is
     the derivative of the biquadratic map through the 9 nodes (xyz) of the
     cell's 2x2 block, at the block coordinates X = (p + xi)/2,
     Y = (q + eta)/2 -- the order 2 analogue of B's bilinear-corner metric
  E  exact velocity projected onto the exact dual basis (control: time and
     hand-over errors only)

Fields must provide basis() / coeffs(t) (trajectory_fields.py) so that all
edge data are precomputed once per time-independent basis field.
"""
import sys
from pathlib import Path

import numpy

sys.path.insert(0, str(Path(__file__).absolute().parent))
from nodal_star import CubedSphere, nodalVectors, greatCircleIntegral, facePoint  # noqa: E402
from cs_chart import _FRAMES, basis as chartBasis  # noqa: E402

ARMS = ('A', 'B', 'C', 'D', 'W', 'R', 'S', 'E')
SNAP = 1e-9


def _geom(f, al, be):
    """Point x and covariant chart basis (a_alpha, a_beta) for per-particle faces f."""
    n, u, w = _FRAMES[f, 0], _FRAMES[f, 1], _FRAMES[f, 2]
    X, Y = numpy.tan(al), numpy.tan(be)
    s = n + X[:, None] * u + Y[:, None] * w
    r = numpy.linalg.norm(s, axis=-1, keepdims=True)
    x = s / r

    def d(ds):
        return (ds - numpy.sum(x * ds, -1, keepdims=True) * x) / r
    return x, d(u / numpy.cos(al)[:, None] ** 2), d(w / numpy.cos(be)[:, None] ** 2)


def _dual(jx, jy):
    g11 = numpy.sum(jx * jx, -1)
    g12 = numpy.sum(jx * jy, -1)
    g22 = numpy.sum(jy * jy, -1)
    det = g11 * g22 - g12 ** 2
    return ((g22[:, None] * jx - g12[:, None] * jy) / det[:, None],
            (-g12[:, None] * jx + g11[:, None] * jy) / det[:, None])


def mintCellIds(M, f, i, j):
    """mint cell Id of cubed-sphere cell (f, i, j), see mintCubedSphereElement."""
    return (f * M + j) * M + i


def mintCubedSphereElement(cs, order):
    """
    mint grid with the topology of the cubed sphere cs, cell (f, i, j) ->
    mintCellIds(M, f, i, j) with its local frame along the panel chart, and
    mint's edge element of the given order on it. The order 2 element uses the
    piecewise cell coordinates of its 2x2 blocks (biquadratic=False), its exact
    geometry in the panel chart where the equiangular grid is uniform.
    """
    import mint
    nid = cs.nodeId
    f2n = numpy.stack([nid[:, :-1, :-1], nid[:, 1:, :-1], nid[:, 1:, 1:], nid[:, :-1, 1:]], -1)
    f2n = numpy.ascontiguousarray(f2n.transpose(0, 2, 1, 3).reshape(-1, 4), numpy.uint64)
    x = cs.nodes
    xyz = numpy.zeros_like(x)
    xyz[:, 0] = numpy.degrees(numpy.arctan2(x[:, 1], x[:, 0]))
    xyz[:, 1] = numpy.degrees(numpy.arcsin(numpy.clip(x[:, 2], -1., 1.)))
    grid = mint.Grid()
    grid.setFlags(1, 1)
    grid.loadFromUgrid2DData(numpy.ascontiguousarray(xyz),
                             f2n, numpy.ascontiguousarray(cs.edges, numpy.uint64))
    # the element keeps the grid alive
    return mint.WhitneyElement(order, grid, biquadratic=False)


def mintCellByCellData(Ex, Ey):
    """
    Cell by cell edge data for mintCubedSphereElement from the circulations
    Ex (..., 6, M, M+1) and Ey (..., 6, M+1, M), oriented +xi, +eta like mint's
    cell edges 0 (south), 1 (east), 2 (north), 3 (west). Returns (..., 6*M*M*4).
    """
    Ex, Ey = numpy.asarray(Ex), numpy.asarray(Ey)
    cbc = numpy.stack([Ex[..., :, :-1], Ey[..., 1:, :], Ex[..., :, 1:], Ey[..., :-1, :]], -1)
    cbc = numpy.swapaxes(cbc, -3, -2)                  # [f, i, j, edge] -> [f, j, i, edge]
    return numpy.ascontiguousarray(cbc.reshape(cbc.shape[:-4] + (-1,)))


def biquadraticBlockBasis(P, f, i, j, xi, eta):
    """
    Covariant basis (dr/dxi, dr/deta) of the biquadratic map through the 9
    nodes P (6, M+1, M+1, 3) of the 2x2 block of cell (f, i, j), at the cell
    coordinates (xi, eta). mint tiles each panel into the blocks of cells
    (2a, 2a + 1) x (2b, 2b + 1), so cell (i, j) is at position
    (p, q) = (i % 2, j % 2) of its block.
    """
    p, q = i % 2, j % 2
    X, Y = 0.5 * (p + xi), 0.5 * (q + eta)

    def lag(t):
        return (numpy.stack([(1 - t) * (1 - 2 * t), 4 * t * (1 - t), t * (2 * t - 1)]),
                numpy.stack([4 * t - 3, 4 - 8 * t, 4 * t - 1]))
    lx, dlx = lag(X)
    ly, dly = lag(Y)
    i0, j0 = i - p, j - q
    dGdX = numpy.zeros(numpy.shape(f) + (3,))
    dGdY = numpy.zeros(numpy.shape(f) + (3,))
    for a in range(3):
        for b in range(3):
            node = P[f, i0 + a, j0 + b]
            dGdX += (dlx[a] * ly[b])[..., None] * node
            dGdY += (lx[a] * dly[b])[..., None] * node
    # dX/dxi = dY/deta = 1/2
    return 0.5 * dGdX, 0.5 * dGdY


class Tracer:

    def __init__(self, M, fld, arms=ARMS, nGauss=12):
        self.M = M
        self.fld = fld
        cs = CubedSphere(M)
        self.cs = cs
        self.a = cs.alphas
        self.da = numpy.diff(cs.alphas)
        P = cs.nodes[cs.nodeId]                          # (6, M+1, M+1, 3)
        self.P = P
        self.normals = _FRAMES[:, 0, :]
        Ex, Ey, Fx, Fy, V, Sx, Sy = [], [], [], [], [], [], []
        for psi, gpsi in fld.basis():
            def vel(x, gpsi=gpsi):
                return numpy.cross(x, gpsi(x))
            Ex.append(greatCircleIntegral(vel, P[:, :-1, :], P[:, 1:, :], nGauss))   # (6, M, M+1)
            Ey.append(greatCircleIntegral(vel, P[:, :, :-1], P[:, :, 1:], nGauss))   # (6, M+1, M)
            if 'W' in arms:
                sx, sy = self._edgeSlopes(vel, nGauss)
                Sx.append(sx)
                Sy.append(sy)
            ps = psi(P)
            Fx.append(ps[:, :, :-1] - ps[:, :, 1:])      # +xi flux through eta-edges (6, M+1, M)
            Fy.append(ps[:, 1:, :] - ps[:, :-1, :])      # +eta flux through xi-edges (6, M, M+1)
            if 'D' in arms:
                Vn = nodalVectors(cs, cs.edgeIntegrals(vel, nGauss))
                V.append(Vn[cs.nodeId])                  # (6, M+1, M+1, 3)
        self.Ex, self.Ey = numpy.array(Ex), numpy.array(Ey)
        self.Sx = numpy.array(Sx) if Sx else None
        self.Sy = numpy.array(Sy) if Sy else None
        self.Fx, self.Fy = numpy.array(Fx), numpy.array(Fy)
        self.V = numpy.array(V) if V else None
        if 'R' in arms or 'S' in arms:
            self._buildMintElement(2)
        self.nSlide = 0
        self.nPanelCross = 0

    def _edgeSlopes(self, vel, nGauss):
        """
        k=2 slope DOFs b = 3 int_0^1 t(s) (2s - 1) ds of the covariant edge
        profile t(s) = v . dx/ds, with s the CELL coordinate (xi along
        xi-edges, eta along eta-edges) of the panel chart -- the same
        coordinate the reconstruction is polynomial in. Oriented like Ex/Ey
        (+xi, +eta); as in whitney2.py, no sign flip is applied to slopes.
        Returns Sx (6, M, M+1), Sy (6, M+1, M).
        """
        sg, wg = numpy.polynomial.legendre.leggauss(nGauss)
        sg, wg = 0.5 * (sg + 1.), 0.5 * wg
        a, da, M = self.a, self.da, self.M
        Sx = numpy.zeros((6, M, M + 1))
        Sy = numpy.zeros((6, M + 1, M))
        for f in range(6):
            for sk, wk in zip(sg, wg):
                # xi-edges: alpha = a_i + s da_i, beta = a_j
                A, B = numpy.meshgrid(a[:-1] + sk * da, a, indexing='ij')
                ja, jb = chartBasis(f, A, B)
                x = facePoint(f, A, B)
                t = numpy.sum(vel(x) * ja, -1) * da[:, None]
                Sx[f] += 3. * wk * t * (2. * sk - 1.)
                # eta-edges: alpha = a_i, beta = a_j + s da_j
                A, B = numpy.meshgrid(a, a[:-1] + sk * da, indexing='ij')
                ja, jb = chartBasis(f, A, B)
                x = facePoint(f, A, B)
                t = numpy.sum(vel(x) * jb, -1) * da[None, :]
                Sy[f] += 3. * wk * t * (2. * sk - 1.)
        return Sx, Sy

    def _buildMintElement(self, order):
        """mint's edge element on the cubed sphere and the cell by cell data of each basis field."""
        self.mintElement = mintCubedSphereElement(self.cs, order)
        self.mintData = mintCellByCellData(self.Ex, self.Ey)

    def _mintCovariant(self, st):
        """Covariant components (v_xi, v_eta) of each basis field from mint's element."""
        f, i, j, xi, eta = st
        cellIds = mintCellIds(self.M, f, i, j)
        ab = [self.mintElement.getCovariantComponents(cellIds, xi, eta, d) for d in self.mintData]
        return numpy.array([a for a, b in ab]), numpy.array([b for a, b in ab])

    def _biquadraticBasis(self, st):
        """Covariant basis of the biquadratic map of each particle's 2x2 block."""
        f, i, j, xi, eta = st
        return biquadraticBlockBasis(self.P, f, i, j, xi, eta)

    # ------------------------------------------------------------------
    def position(self, st):
        f, i, j, xi, eta = st
        al = self.a[i] + xi * self.da[i]
        be = self.a[j] + eta * self.da[j]
        return _geom(f, al, be)[0]

    def locate(self, x, exclude=None):
        """xyz -> state; exclude: per-particle face to avoid (panel hand-over)."""
        proj = x @ self.normals.T
        if exclude is not None:
            proj[numpy.arange(len(x)), exclude] = -numpy.inf
        f = numpy.argmax(proj, axis=-1)
        n, u, w = _FRAMES[f, 0], _FRAMES[f, 1], _FRAMES[f, 2]
        xn = numpy.sum(x * n, -1)
        lim = numpy.pi / 4
        al = numpy.clip(numpy.arctan(numpy.sum(x * u, -1) / xn), -lim, lim)
        be = numpy.clip(numpy.arctan(numpy.sum(x * w, -1) / xn), -lim, lim)
        i = numpy.clip(numpy.searchsorted(self.a, al) - 1, 0, self.M - 1)
        j = numpy.clip(numpy.searchsorted(self.a, be) - 1, 0, self.M - 1)
        xi = numpy.clip((al - self.a[i]) / self.da[i], 0., 1.)
        eta = numpy.clip((be - self.a[j]) / self.da[j], 0., 1.)
        for c in (xi, eta):
            c[c < SNAP] = 0.
            c[c > 1 - SNAP] = 1.
        return [f, i, j, xi, eta]

    # ------------------------------------------------------------------
    def velocity(self, arm, st, t):
        """Cell-local contravariant velocity (v^xi, v^eta) at the particles."""
        f, i, j, xi, eta = st
        al = self.a[i] + xi * self.da[i]
        be = self.a[j] + eta * self.da[j]
        x, ja, jb = _geom(f, al, be)
        jx = ja * self.da[i][:, None]
        jy = jb * self.da[j][:, None]
        if arm == 'E':
            v = self.fld.velocity(x, t)
            dx, dy = _dual(jx, jy)
            return numpy.sum(v * dx, -1), numpy.sum(v * dy, -1)
        c = self.fld.coeffs(t)                   # (nbasis,) + t.shape
        if arm in 'ABWRS':
            if arm in 'RS':
                avxi, aveta = self._mintCovariant(st)
                vxi = numpy.sum(c * avxi, 0)
                veta = numpy.sum(c * aveta, 0)
            else:
                vxi = numpy.sum(c * ((1 - eta) * self.Ex[:, f, i, j] + eta * self.Ex[:, f, i, j + 1]), 0)
                veta = numpy.sum(c * ((1 - xi) * self.Ey[:, f, i, j] + xi * self.Ey[:, f, i + 1, j]), 0)
            if arm == 'W':
                vxi = vxi + (2 * xi - 1) * numpy.sum(
                    c * ((1 - eta) * self.Sx[:, f, i, j] + eta * self.Sx[:, f, i, j + 1]), 0)
                veta = veta + (2 * eta - 1) * numpy.sum(
                    c * ((1 - xi) * self.Sy[:, f, i, j] + xi * self.Sy[:, f, i + 1, j]), 0)
            if arm == 'S':
                jx, jy = self._biquadraticBasis(st)
            if arm == 'B':
                P = self.P
                p00, p10, p01, p11 = P[f, i, j], P[f, i + 1, j], P[f, i, j + 1], P[f, i + 1, j + 1]
                jx = (1 - eta)[:, None] * (p10 - p00) + eta[:, None] * (p11 - p01)
                jy = (1 - xi)[:, None] * (p01 - p00) + xi[:, None] * (p11 - p10)
            g11 = numpy.sum(jx * jx, -1)
            g12 = numpy.sum(jx * jy, -1)
            g22 = numpy.sum(jy * jy, -1)
            det = g11 * g22 - g12 ** 2
            return (g22 * vxi - g12 * veta) / det, (-g12 * vxi + g11 * veta) / det
        if arm == 'C':
            J = numpy.sum(x * numpy.cross(jx, jy), -1)
            fxi = numpy.sum(c * ((1 - xi) * self.Fx[:, f, i, j] + xi * self.Fx[:, f, i + 1, j]), 0)
            feta = numpy.sum(c * ((1 - eta) * self.Fy[:, f, i, j] + eta * self.Fy[:, f, i, j + 1]), 0)
            return fxi / J, feta / J
        if arm == 'D':
            V = self.V
            w00 = ((1 - xi) * (1 - eta))[None, :, None]
            w10 = (xi * (1 - eta))[None, :, None]
            w01 = ((1 - xi) * eta)[None, :, None]
            w11 = (xi * eta)[None, :, None]
            vb = (w00 * V[:, f, i, j] + w10 * V[:, f, i + 1, j]
                  + w01 * V[:, f, i, j + 1] + w11 * V[:, f, i + 1, j + 1])
            v = numpy.sum(c[..., None] * vb, 0)
            dx, dy = _dual(jx, jy)
            return numpy.sum(v * dx, -1), numpy.sum(v * dy, -1)
        raise ValueError(arm)

    # ------------------------------------------------------------------
    def _neighbour(self, st, axis, side):
        """Hand particles sitting on cell face (axis, side) over to the adjacent cell."""
        f, i, j, xi, eta = (a.copy() for a in st)
        idx = i if axis == 0 else j
        new = idx + (1 if side == 1 else -1)
        inside = (new >= 0) & (new < self.M)
        if axis == 0:
            i[inside] = new[inside]
            xi[inside] = 1. - side
        else:
            j[inside] = new[inside]
            eta[inside] = 1. - side
        out = ~inside
        if numpy.any(out):
            self.nPanelCross += int(out.sum())
            x = self.position([s[out] for s in st])
            g = self.locate(x, exclude=st[0][out])
            for a, b in zip((f, i, j, xi, eta), g):
                a[out] = b
        return [f, i, j, xi, eta]

    @staticmethod
    def _onFace(st, axis):
        c = st[3 + axis]
        return numpy.where(c == 0., 0, numpy.where(c == 1., 1, -1))

    def _handover(self, arm, st, t):
        """
        At the start of a step: particles on a cell face whose velocity points
        out of the cell move to the neighbour, unless the neighbour's velocity
        points straight back (converging normal velocity, Whitney) -> slide.
        Returns the new state and per-particle slide masks for (xi, eta).
        """
        slide = [numpy.zeros(len(st[0]), bool), numpy.zeros(len(st[0]), bool)]
        for axis in (0, 1):
            side = self._onFace(st, axis)
            if not numpy.any(side >= 0):
                continue
            v = self.velocity(arm, st, t)[axis]
            outw = ((side == 1) & (v > 0)) | ((side == 0) & (v < 0))
            for s in (0, 1):
                m = outw & (side == s)
                if not numpy.any(m):
                    continue
                sub = [a[m] for a in st]
                nb = self._neighbour(sub, axis, s)
                vn = self.velocity(arm, nb, t[m] if numpy.ndim(t) else t)
                # does the neighbour's velocity point straight back across the SAME face?
                # (only the crossed face counts: at a cell corner the particle may
                # legitimately be leaving the neighbour through its other face)
                back = numpy.zeros(m.sum(), bool)
                same = nb[0] == sub[0]
                for ax2 in (0, 1):
                    c2, idx2 = nb[3 + ax2], nb[1 + ax2]
                    if ax2 == axis:
                        crossed = same.copy()
                    else:
                        crossed = numpy.zeros_like(same)
                    # across panels: the crossed face is the one on the panel boundary
                    pb = ((idx2 == 0) & (c2 == 0.)) | ((idx2 == self.M - 1) & (c2 == 1.))
                    crossed |= (~same) & pb
                    sd2 = numpy.where(c2 == 0., 0, 1)
                    back |= crossed & (((sd2 == 1) & (vn[ax2] > 0)) | ((sd2 == 0) & (vn[ax2] < 0)))
                acc = ~back
                mi = numpy.where(m)[0]
                for a, b in zip(st, nb):
                    a[mi[acc]] = b[acc]
                slide[axis][mi[back]] = True
                self.nSlide += int(back.sum())
        return st, slide

    def _rk4(self, arm, st, t, h, slide):
        f, i, j, xi, eta = st

        def vel(xx, ee, tt):
            vx, ve = self.velocity(arm, [f, i, j, xx, ee], tt)
            return numpy.where(slide[0], 0., vx), numpy.where(slide[1], 0., ve)
        k1 = vel(xi, eta, t)
        k2 = vel(xi + 0.5 * h * k1[0], eta + 0.5 * h * k1[1], t + 0.5 * h)
        k3 = vel(xi + 0.5 * h * k2[0], eta + 0.5 * h * k2[1], t + 0.5 * h)
        k4 = vel(xi + h * k3[0], eta + h * k3[1], t + h)
        nx = xi + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        ne = eta + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        return nx, ne, k1, vel(nx, ne, t + h)

    def advance(self, arm, st, t0, t1, dt, visited=None):
        """
        Advance all particles from t0 to t1 (in place); returns st.
        visited: optional bool array, OR-ed with "particle is in a cube-corner cell".
        """
        N = len(st[0])
        t = numpy.full(N, float(t0))
        for _ in range(100 * int(numpy.ceil((t1 - t0) / dt)) + 1000):
            act = numpy.where(t1 - t > 1e-14)[0]
            if len(act) == 0:
                return st
            sub = [a[act] for a in st]
            ts = t[act]
            sub, slide = self._handover(arm, sub, ts)
            h = numpy.minimum(dt, t1 - ts)
            nx, ne, k0, k1 = self._rk4(arm, sub, ts, h, slide)
            # exit detection on the non-sliding axes
            sstar = numpy.ones(len(act))
            for axis, (y0, y1, f0, f1) in enumerate(((sub[3], nx, k0[0], k1[0]),
                                                     (sub[4], ne, k0[1], k1[1]))):
                for b in (0., 1.):
                    viol = (~slide[axis]) & ((y1 < b) if b == 0. else (y1 > b))
                    if not numpy.any(viol):
                        continue
                    s = _hermiteRoot(y0[viol], y1[viol], h[viol] * f0[viol], h[viol] * f1[viol], b)
                    sstar[viol] = numpy.minimum(sstar[viol], s)
            ex = sstar < 1.
            if numpy.any(ex):
                hx = h[ex] * sstar[ex]
                sx = [a[ex] for a in sub]
                slx = [slide[0][ex], slide[1][ex]]
                ex_x, ex_e, _, _ = self._rk4(arm, sx, ts[ex], hx, slx)
                nx[ex] = ex_x
                ne[ex] = ex_e
                h[ex] = hx
            for c in (nx, ne):
                c[:] = numpy.clip(c, 0., 1.)
                c[ex & (c < SNAP)] = 0.
                c[ex & (c > 1 - SNAP)] = 1.
            sub[3], sub[4] = nx, ne
            for a, b in zip(st, sub):
                a[act] = b
            t[act] = ts + h
            if visited is not None:
                edge = (0, self.M - 1)
                visited[act] |= numpy.isin(sub[1], edge) & numpy.isin(sub[2], edge)
        raise RuntimeError('advance: too many iterations')


def _hermiteRoot(y0, y1, d0, d1, b, nSample=32, nBisect=50):
    """
    First s in (0, 1] at which the cubic Hermite H(s) (H(0)=y0, H(1)=y1,
    H'=d0, d1) leaves the cell through the face at b. The cubic need not be
    monotone (e.g. a particle on or near the face moving inwards first), so
    the first crossing is bracketed by sampling and then bisected.
    """
    sgn = 1. if b == 1. else -1.

    def out(s):
        s2, s3 = s * s, s * s * s
        H = ((2 * s3 - 3 * s2 + 1) * y0 + (s3 - 2 * s2 + s) * d0
             + (-2 * s3 + 3 * s2) * y1 + (s3 - s2) * d1)
        return sgn * (H - b) > 0
    sg = numpy.linspace(0., 1., nSample + 1)
    outs = numpy.stack([out(numpy.full_like(y0, sk)) for sk in sg[1:]])   # (nSample, n)
    outs[-1] = True                                   # y1 is outside by construction
    k = numpy.argmax(outs, axis=0)                     # first sample outside
    lo, hi = sg[k], sg[k + 1]
    for _ in range(nBisect):
        mid = 0.5 * (lo + hi)
        o = out(mid)
        hi = numpy.where(o, mid, hi)
        lo = numpy.where(o, lo, mid)
    return hi
