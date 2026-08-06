#!/usr/bin/env python
"""
Convergence of MINT's line-integral error for the pole-asymmetric stream
function

    psi(lambda, theta) = cos(theta) * (1 + sin(theta)) * cos(lambda)

on latitude-longitude AND cubed-sphere grids, at increasing resolution, for
several test lines with varying proximity to the poles -- see
mint/tests/test_pole_asymmetric.py for the analytic background. This is the
azimuthal wavenumber-1 mode with envelope g(sin(theta)) = 1 + sin(theta),
giving finite nonzero velocity at the north pole (g(1) = 2) and exactly
zero velocity at the south pole (g(-1) = 0).

For each (grid family, test line) combination, the flux computed by
mint.PolylineIntegral is compared against the EXACT line integral (not a
quadrature approximation -- see exactLineIntegral below) at increasing grid
resolution, and |error| is plotted against N = number of grid cells on
log-log axes, alongside a reference N^-1 slope.

Usage: python scripts/pole_asymmetric_convergence.py
"""
import numpy
import sympy
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mint import Grid, PolylineIntegral, NUM_EDGES_PER_QUAD, CELL_BY_CELL_DATA

HERE = Path(__file__).absolute().parent
DATA_DIR = HERE.parent / 'data'
DEG2RAD = numpy.pi / 180.0


# ---------------------------------------------------------------------------
# The pole-asymmetric field (see mint/tests/test_pole_asymmetric.py)
# ---------------------------------------------------------------------------
def uDotDl(lam, theta, dlam, dtheta):
    """u . dl integrand at (lam, theta) for a local displacement (dlam, dtheta)."""
    return (numpy.cos(theta) * (numpy.cos(2. * theta) - numpy.sin(theta)) * numpy.cos(lam) * dlam
            + (1. + numpy.sin(theta)) * numpy.sin(lam) * dtheta)


def straightLineIntegral(p0, p1, n=50):
    """
    Line integral of u . dl along the path STRAIGHT IN (lon, lat) from p0
    to p1 (lon, lat in degrees; trailing axis of size 2, batchable),
    evaluated by an n-point quadrature. Used only for buildEdgeData below,
    where it is called thousands of times on short grid edges (quadrature
    error there is utterly negligible relative to the errors being
    studied); see exactLineIntegral for the reference value used to judge
    those errors, which has no quadrature error at all.
    """
    p0 = numpy.asarray(p0, dtype=numpy.float64)
    p1 = numpy.asarray(p1, dtype=numpy.float64)
    lam0, th0 = p0[..., 0] * DEG2RAD, p0[..., 1] * DEG2RAD
    lam1, th1 = p1[..., 0] * DEG2RAD, p1[..., 1] * DEG2RAD
    dlam, dth = lam1 - lam0, th1 - th0

    tt = numpy.linspace(0., 1., n + 1)
    lam = lam0[..., None] + tt * dlam[..., None]
    theta = th0[..., None] + tt * dth[..., None]
    integrand = uDotDl(lam, theta, dlam[..., None], dth[..., None])

    dt = tt[1] - tt[0]
    return dt * (integrand[..., 0] / 2. + integrand[..., 1:-1].sum(axis=-1) + integrand[..., -1] / 2.)


# ---------------------------------------------------------------------------
# EXACT (not quadrature) straight-line integral, used only for the handful
# of whole-polyline reference values (one per test line -- cheap even
# though sympy is slow per call). Substituting the straight line into
# u . dl turns it into a finite trigonometric polynomial in t (products and
# sums of sin/cos of expressions linear in t); rewriting sin/cos as complex
# exponentials turns every term into coeff * exp(intercept + slope * t),
# each of which integrates over [0, 1] in closed form -- see
# work/vector_potential/vector_potential.py's
# _integrate_unit_interval_trig_polynomial for the derivation and further
# discussion of why this is both exact and much faster than sympy's
# general-purpose integrator on an expression like this.
# ---------------------------------------------------------------------------
_lam_s, _th_s, _t_s, _dlam_s, _dth_s = sympy.symbols('lam theta t dlam dtheta', real=True)
_UDOTDL_SYM = (sympy.cos(_th_s) * (sympy.cos(2 * _th_s) - sympy.sin(_th_s)) * sympy.cos(_lam_s) * _dlam_s
               + (1 + sympy.sin(_th_s)) * sympy.sin(_lam_s) * _dth_s)


def _integrateUnitIntervalTrigPolynomial(expr, t_sym):
    expr_exp = sympy.expand(expr.rewrite(sympy.exp))
    total = sympy.Integer(0)
    for term in sympy.Add.make_args(expr_exp):
        factors = term.as_ordered_factors()
        exp_args = [f.args[0] for f in factors if f.func == sympy.exp]
        rest = [f for f in factors if f.func != sympy.exp]
        total_arg = sympy.expand(sum(exp_args)) if exp_args else sympy.Integer(0)
        coeff = sympy.Mul(*rest) if rest else sympy.Integer(1)
        slope = sympy.diff(total_arg, t_sym)
        intercept = total_arg.subs(t_sym, 0)
        if slope == 0:
            total += coeff * sympy.exp(intercept)
        else:
            total += coeff * sympy.exp(intercept) * (sympy.exp(slope) - 1) / slope
    return sympy.re(sympy.expand_complex(total))


def exactLineIntegral(p0, p1):
    """
    EXACT (no quadrature) line integral of u . dl along the path straight
    in (lon, lat) from p0 to p1 (lon, lat in degrees).
    """
    lamA, thA = p0[0] * DEG2RAD, p0[1] * DEG2RAD
    lamB, thB = p1[0] * DEG2RAD, p1[1] * DEG2RAD
    dlam, dth = lamB - lamA, thB - thA

    lam_t = lamA + _t_s * dlam
    th_t = thA + _t_s * dth
    integrand = sympy.expand(_UDOTDL_SYM.subs({_lam_s: lam_t, _th_s: th_t, _dlam_s: dlam, _dth_s: dth}))
    return float(_integrateUnitIntervalTrigPolynomial(integrand, _t_s))


def buildEdgeData(grid):
    """
    Quadrature-exact cell-by-cell edge-integrated data for u on the grid's
    own edges. Note: this assumes each edge is straight in (lon, lat),
    which is exact for the lat-lon grids below but only approximate for
    the cubed-sphere grids (whose true edges are closer to great-circle
    arcs) -- that extra approximation is itself a discretization error
    that should shrink with resolution, same as the rest of the scheme.
    """
    ncells = grid.getNumberOfCells()
    points = grid.getPoints()  # (ncells, 4, 3): lon, lat, elev in degrees
    data = numpy.zeros((ncells, NUM_EDGES_PER_QUAD))
    for i0 in range(NUM_EDGES_PER_QUAD):
        i1 = (i0 + 1) % NUM_EDGES_PER_QUAD
        # edges 2 and 3 run backwards relative to the (i0 -> i1) node order
        sign = 1 - 2 * (i0 // 2)
        data[:, i0] = sign * straightLineIntegral(points[:, i0, :2], points[:, i1, :2])
    return data


# ---------------------------------------------------------------------------
# Grid families. Both are loaded with coordinates in degrees (confirmed by
# inspection: the "cs_N.nc$physics" cubed-sphere files give lon in
# [0, 360], lat in [-90, 90], unlike e.g. "cubedsphereN.nc$grid" or
# "mesh_CN.nc$unit_test", which are in radians regardless of the `degrees`
# flag passed to Grid.setFlags).
# ---------------------------------------------------------------------------
def loadLatLon(n):
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=True)
    grid.loadFromUgrid2DFile(f'{DATA_DIR}/latlon{n}x{n // 2}Shifted.nc$mesh')
    return grid


def loadCubedSphere(n):
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=1, averageLonAtPole=1, degrees=True)
    grid.loadFromUgrid2DFile(f'{DATA_DIR}/cs_{n}.nc$physics')
    return grid


# (name, resolutions, loader, periodX)
GRID_FAMILIES = [
    ('lon-lat', [20, 40, 80, 160], loadLatLon, 360.),
    ('cubed-sphere', [4, 16, 64], loadCubedSphere, 360.),
]


# ---------------------------------------------------------------------------
# Test lines: (name, p0, p1), all straight in (lon, lat). Some end at or
# very close to a pole, on purpose.
# ---------------------------------------------------------------------------
TEST_LINES = [
    ('away from poles',          (-30.0, -20.0), (45.0, 50.0)),
    ('1 deg from north pole',    (-30.0, 60.0),  (45.0, 89.0)),
    ('AT north pole (finite u)', (-30.0, -20.0), (45.0, 90.0)),
    ('AT south pole (zero u)',   (-30.0, -20.0), (45.0, -90.0)),
]


def computeFlux(grid, periodX, p0, p1, data):
    xyz = numpy.array([(p0[0], p0[1], 0.), (p1[0], p1[1], 0.)])
    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=128, periodX=periodX, enableFolding=False)
    pli.computeWeights(xyz, counterclock=False)
    return pli.getIntegral(data, placement=CELL_BY_CELL_DATA)


def main():
    # results[(familyName, lineName)] = (list of numCells, list of |error|)
    results = {}

    for familyName, resolutions, loader, periodX in GRID_FAMILIES:
        for lineName, p0, p1 in TEST_LINES:
            exact = exactLineIntegral(p0, p1)
            numCellsList, errs = [], []
            for n in resolutions:
                grid = loader(n)
                numCells = grid.getNumberOfCells()
                data = buildEdgeData(grid)
                flux = computeFlux(grid, periodX, p0, p1, data)
                err = abs(flux - exact)
                numCellsList.append(numCells)
                errs.append(err)
                print(f'{familyName:13s} | {lineName:26s} | ncells={numCells:7d} '
                      f'flux={flux:12.6f} exact={exact:12.6f} err={err:.3e}')
            results[(familyName, lineName)] = (numCellsList, errs)
        print()

    plotResults(results)


def plotResults(results):
    familyNames = sorted({k[0] for k in results})
    fig, axes = plt.subplots(1, len(familyNames), figsize=(7 * len(familyNames), 6), squeeze=False)
    axes = axes[0]

    colors = plt.cm.tab10.colors

    for ax, familyName in zip(axes, familyNames):
        for i, (lineName, _, _) in enumerate(TEST_LINES):
            numCellsList, errs = results[(familyName, lineName)]
            ax.loglog(numCellsList, errs, 'o-', color=colors[i % len(colors)], label=lineName)

        # reference N^-1 slope, anchored on the first test line's first point
        anchorN, anchorErr = results[(familyName, TEST_LINES[0][0])][0][0], \
            results[(familyName, TEST_LINES[0][0])][1][0]
        allN = sorted({n for lineName, _, _ in TEST_LINES for n in results[(familyName, lineName)][0]})
        refN = numpy.array([allN[0], allN[-1]], dtype=numpy.float64)
        refErr = anchorErr * (refN / anchorN) ** (-1.0)
        ax.loglog(refN, refErr, 'k--', linewidth=1.5, label=r'$N^{-1}$ reference')

        ax.set_xlabel('N = number of grid cells')
        ax.set_ylabel('|line integral error|')
        ax.set_title(familyName)
        ax.legend(fontsize='small')
        ax.grid(True, which='both', alpha=0.3)

    fig.suptitle(r'Convergence of $\int u\cdot dl$ for $\psi=\cos\theta(1+\sin\theta)\cos\lambda$')
    fig.tight_layout()
    outfile = HERE / 'pole_asymmetric_convergence.png'
    fig.savefig(outfile, dpi=150)
    print(f'wrote {outfile}')


if __name__ == '__main__':
    main()
