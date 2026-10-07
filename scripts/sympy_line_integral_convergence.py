#!/usr/bin/env python
"""
Symbolic estimate of the convergence rate of a straight-line integral of a
1-form reconstructed with the lowest-order mimetic (Whitney-type) edge basis
functions on a square cell of size h x h centred at (0, 0).

Cell parametric coordinates 0 <= xi, eta <= 1:   x = h (xi - 1/2),  y = h (eta - 1/2)

Basis 1-forms (in parametric coordinates):
    bottom : (1 - eta) d xi
    east   : xi d eta
    north  : eta d xi
    west   : (1 - xi) d eta
(edges oriented in the +xi / +eta directions).

Steps
  1. interpolation coefficients = line integral of the 1-form over each edge
  2. straight target line from (xi0, eta0) to (xi1, eta1)
  3. approximate integral = sum_e coeff_e * int_line basis_e
  4. compare with the exact line integral of the (same) truncated 1-form;
     the leading power of h in the difference is the convergence rate.

The 1-form  w = u dx + v dy  is Taylor expanded about (0,0) up to and including
quadratic terms in x, y (cubic and higher terms are dropped).
"""
import sympy as sp

h = sp.symbols('h', positive=True)
xi, eta, t = sp.symbols('xi eta t', real=True)
x, y = sp.symbols('x y', real=True)

# ---- Taylor coefficients of the 1-form w = u dx + v dy -----------------------
u0, ux, uy, uxx, uxy, uyy = sp.symbols('u0 u_x u_y u_xx u_xy u_yy', real=True)
v0, vx, vy, vxx, vxy, vyy = sp.symbols('v0 v_x v_y v_xx v_xy v_yy', real=True)


def taylor2(c0, cx, cy, cxx, cxy, cyy):
    return (c0 + cx*x + cy*y
            + sp.Rational(1, 2)*cxx*x**2 + cxy*x*y + sp.Rational(1, 2)*cyy*y**2)


u = taylor2(u0, ux, uy, uxx, uxy, uyy)
v = taylor2(v0, vx, vy, vxx, vxy, vyy)

# parametric -> physical
xmap = h*(xi - sp.Rational(1, 2))
ymap = h*(eta - sp.Rational(1, 2))
dxdxi, dydeta = h, h

# pull-back of w to parametric coordinates: w = a d xi + b d eta
a = sp.expand(u.subs({x: xmap, y: ymap}, simultaneous=True)*dxdxi)
b = sp.expand(v.subs({x: xmap, y: ymap}, simultaneous=True)*dydeta)


def line_integral(a_expr, b_expr, p0, p1):
    """int of a d xi + b d eta along the straight line p0 -> p1."""
    sub = {xi: p0[0] + t*(p1[0] - p0[0]), eta: p0[1] + t*(p1[1] - p0[1])}
    dxi, deta = p1[0] - p0[0], p1[1] - p0[1]
    integrand = a_expr.subs(sub, simultaneous=True)*dxi \
        + b_expr.subs(sub, simultaneous=True)*deta
    return sp.integrate(sp.expand(integrand), (t, 0, 1))


# ---- Step 1: interpolation coefficients (edge line integrals) ---------------
c_bottom = line_integral(a, b, (0, 0), (1, 0))
c_east = line_integral(a, b, (1, 0), (1, 1))
c_north = line_integral(a, b, (0, 1), (1, 1))
c_west = line_integral(a, b, (0, 0), (0, 1))

# basis 1-forms as (coef of d xi, coef of d eta)
basis = {
    'bottom': ((1 - eta), 0, c_bottom),
    'east':   (0, xi, c_east),
    'north':  (eta, 0, c_north),
    'west':   (0, (1 - xi), c_west),
}


def convergence(p0, p1, verbose=True):
    """Return (approx, exact, error, leading order in h of the error)."""
    # Step 3: sum of coefficients times line integral of each basis function
    approx = 0
    for name, (ba, bb, c) in basis.items():
        approx += c*line_integral(sp.sympify(ba), sp.sympify(bb), p0, p1)
    approx = sp.expand(approx)
    # Step 4: exact line integral of the truncated 1-form
    exact = sp.expand(line_integral(a, b, p0, p1))
    err = sp.expand(approx - exact)
    if err == 0:
        return approx, exact, err, sp.oo
    poly = sp.Poly(err, h)
    order = min(m[0] for m in poly.monoms())
    return approx, exact, err, order


def report(label, p0, p1):
    approx, exact, err, order = convergence(p0, p1)
    print(f'--- {label}: ({p0}) -> ({p1})')
    if err == 0:
        print('    error is identically zero (exact for quadratic 1-forms)')
        return
    lead = sp.factor(sp.Poly(err, h).coeff_monomial(h**order))
    print(f'    error = O(h^{order});  exact integral = O(h^'
          f'{min(m[0] for m in sp.Poly(exact, h).monoms())})')
    print(f'    leading error coeff (h^{order}): {lead}')
    # remaining powers
    print(f'    full error: {sp.collect(err, h)}')


if __name__ == '__main__':
    print('Interpolation coefficients (edge line integrals):')
    for n, c in [('bottom', c_bottom), ('east', c_east),
                 ('north', c_north), ('west', c_west)]:
        print(f'  {n:6s}: {sp.collect(sp.expand(c), h)}')
    print()

    R = sp.Rational
    # generic line with symbolic endpoints
    s0, s1, s2, s3 = sp.symbols('xi0 eta0 xi1 eta1', real=True)
    approx, exact, err, order = convergence((s0, s1), (s2, s3))
    print('Generic line (xi0,eta0)->(xi1,eta1):')
    print(f'  error = O(h^{order})')
    print(f'  leading coefficient:')
    print('   ', sp.factor(sp.Poly(err, h).coeff_monomial(h**order)))
    print()

    lines = {
        'bottom edge':        ((0, 0), (1, 0)),
        'east edge':          ((1, 0), (1, 1)),
        'horizontal midline': ((0, R(1, 2)), (1, R(1, 2))),
        'vertical midline':   ((R(1, 2), 0), (R(1, 2), 1)),
        'diagonal':           ((0, 0), (1, 1)),
        'anti-diagonal':      ((0, 1), (1, 0)),
        'off-centre line':    ((R(1, 5), R(1, 10)), (R(9, 10), R(7, 10))),
        'short interior':     ((R(1, 4), R(1, 4)), (R(3, 4), R(1, 2))),
    }
    for label, (p0, p1) in lines.items():
        report(label, p0, p1)
