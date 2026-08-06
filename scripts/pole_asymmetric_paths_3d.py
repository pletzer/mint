#!/usr/bin/env python
"""
3-D visualization (pyvista) of the pole-asymmetric velocity field

    psi(lambda, theta) = cos(theta) * (1 + sin(theta)) * cos(lambda)

on the sphere, together with the straight-in-(lon, lat) integration paths
used by pole_asymmetric_convergence.py and mint/tests/test_pole_asymmetric.py
(TEST_LINES is imported from pole_asymmetric_convergence.py so both scripts
always agree on exactly which paths are being studied).

Arrows are drawn at their TRUE magnitude (not normalized), colouring the
sphere surface by speed too -- as established in the accompanying
work/vector_potential visualizations, this is the only way to see the pole
regularity property visually: the north pole has a finite, nonzero limit
(g(1) = 2) while the south pole is exactly zero (g(-1) = 0), and each
straight-in-(lon, lat) path is drawn as a 3-D curve on the sphere (it is
only straight in the (lon, lat) chart, not in 3-D) with its two endpoints
marked, so the paths ending at/near each pole can be seen against the field
they are integrating.

Usage: python scripts/pole_asymmetric_paths_3d.py
"""
import sys
import numpy
import pyvista as pv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parent))
from pole_asymmetric_convergence import TEST_LINES, DEG2RAD  # noqa: E402

R = 1.0


def velocityCartesian(lam, theta):
    """
    Tangent-plane velocity of the pole-asymmetric field, in Cartesian
    components. Returns (X, Y, Z, U, V, W): (X, Y, Z) points on the unit
    sphere, (U, V, W) the velocity vector there.
    """
    lam = numpy.asarray(lam, dtype=numpy.float64)
    theta = numpy.asarray(theta, dtype=numpy.float64)

    u_lam = (numpy.cos(2. * theta) - numpy.sin(theta)) * numpy.cos(lam)
    u_th = (1. + numpy.sin(theta)) * numpy.sin(lam)

    X = R * numpy.cos(theta) * numpy.cos(lam)
    Y = R * numpy.cos(theta) * numpy.sin(lam)
    Z = R * numpy.sin(theta)

    e_lam = numpy.stack([-numpy.sin(lam), numpy.cos(lam), numpy.zeros_like(lam)], axis=-1)
    e_th = numpy.stack([-numpy.sin(theta) * numpy.cos(lam),
                         -numpy.sin(theta) * numpy.sin(lam),
                         numpy.cos(theta)], axis=-1)

    vel = u_lam[..., None] * e_lam + u_th[..., None] * e_th
    return X, Y, Z, vel[..., 0], vel[..., 1], vel[..., 2]


def pathToCartesian(p0, p1, n=200):
    """Sample the straight-in-(lon, lat) path from p0 to p1 (degrees) at
    n+1 points and return its 3-D Cartesian coordinates, shape (n+1, 3)."""
    t = numpy.linspace(0., 1., n + 1)
    lam = (p0[0] + t * (p1[0] - p0[0])) * DEG2RAD
    theta = (p0[1] + t * (p1[1] - p0[1])) * DEG2RAD
    x = R * numpy.cos(theta) * numpy.cos(lam)
    y = R * numpy.cos(theta) * numpy.sin(lam)
    z = R * numpy.sin(theta)
    return numpy.column_stack([x, y, z])


def buildArrowPoints():
    """Coarse global grid plus dense near-pole rings, matching the style
    of work/vector_potential/pyvista_pole_asymmetric.py."""
    points = []
    n_lat, n_lon = 31, 32
    lat_vals = numpy.linspace(-numpy.pi / 2, numpy.pi / 2, n_lat)[1:-1]
    lon_vals = numpy.linspace(-numpy.pi, numpy.pi, n_lon, endpoint=False)
    for th in lat_vals:
        for lm in lon_vals:
            points.append((lm, th))

    n_ring = 48
    lon_ring = numpy.linspace(-numpy.pi, numpy.pi, n_ring, endpoint=False)
    for eps_deg in (8, 4, 2, 1, 0.5, 0.25, 0.1):
        eps = numpy.radians(eps_deg)
        for sign in (+1, -1):
            th = sign * (numpy.pi / 2 - eps)
            for lm in lon_ring:
                points.append((lm, th))

    lam, theta = zip(*points)
    return numpy.array(lam), numpy.array(theta)


def main():
    lam, theta = buildArrowPoints()
    X, Y, Z, U, V, W = velocityCartesian(lam, theta)
    speed = numpy.sqrt(U**2 + V**2 + W**2)

    cloud = pv.PolyData(numpy.column_stack([X, Y, Z]))
    cloud['vectors'] = numpy.column_stack([U, V, W])
    cloud['speed'] = speed
    # true magnitude, global scale (arrows vanish near the south pole, stay
    # finite near the north pole -- see work/vector_potential's writeup)
    factor = 0.3 / speed.max()
    glyphs = cloud.glyph(orient='vectors', scale='speed', factor=factor)

    sphere_mesh = pv.Sphere(radius=R, theta_resolution=180, phi_resolution=360)
    pts = sphere_mesh.points
    lam_m = numpy.arctan2(pts[:, 1], pts[:, 0])
    theta_m = numpy.arcsin(numpy.clip(pts[:, 2] / R, -1.0, 1.0))
    _, _, _, Um, Vm, Wm = velocityCartesian(lam_m, theta_m)
    sphere_mesh['speed'] = numpy.sqrt(Um**2 + Vm**2 + Wm**2)

    # 'white' would be invisible against the white legend background used
    # below, hence 'black' instead
    path_colors = ['black', 'orange', 'magenta', 'cyan']

    def new_plotter(title, poles='both'):
        pl = pv.Plotter(off_screen=True, window_size=(1300, 1300))
        pl.add_mesh(sphere_mesh, scalars='speed', cmap='viridis',
                    show_scalar_bar=True, scalar_bar_args={'title': 'speed |u|'})
        pl.add_mesh(glyphs, scalars='speed', cmap='viridis', show_scalar_bar=False)

        for (lineName, p0, p1), color in zip(TEST_LINES, path_colors):
            curve = pathToCartesian(p0, p1)
            pl.add_mesh(pv.MultipleLines(points=curve).tube(radius=0.006), color=color, label=lineName)
            pl.add_mesh(pv.Sphere(radius=0.015, center=curve[0]), color=color)
            pl.add_mesh(pv.Sphere(radius=0.015, center=curve[-1]), color=color)

        # only label the pole(s) relevant to this view -- with both always
        # on, a close-up centred on one pole would still show the other
        # pole's label projected somewhere in the same frame, which reads
        # as a mislabeled marker
        if poles in ('both', 'north'):
            pl.add_mesh(pv.Sphere(radius=0.02 * R, center=(0, 0, R)), color='red')
            pl.add_point_labels([[0, 0, 1.15 * R]], ['N'], font_size=24, text_color='red',
                                  shape=None, always_visible=True)
        if poles in ('both', 'south'):
            pl.add_mesh(pv.Sphere(radius=0.02 * R, center=(0, 0, -R)), color='blue')
            pl.add_point_labels([[0, 0, -1.15 * R]], ['S'], font_size=24, text_color='blue',
                                  shape=None, always_visible=True)

        pl.add_legend(bcolor='white', face=None)
        pl.add_text(title, font_size=11)
        return pl

    pl = new_plotter('psi = cos(theta)(1+sin(theta))cos(lambda): velocity field and integration paths')
    pl.camera_position = [(2.6, 2.6, 1.8), (0, 0, 0), (0, 0, 1)]
    pl.screenshot(str(Path(__file__).parent / 'pole_asymmetric_paths_3d_overview.png'))
    print('wrote pole_asymmetric_paths_3d_overview.png')

    pl_north = new_plotter('North pole close-up: path ending here has finite velocity', poles='north')
    pl_north.camera_position = [(0, 0.01, 3.2), (0, 0, 1), (0, 1, 0)]
    pl_north.screenshot(str(Path(__file__).parent / 'pole_asymmetric_paths_3d_north_pole.png'))
    print('wrote pole_asymmetric_paths_3d_north_pole.png')

    pl_south = new_plotter('South pole close-up: path ending here has zero velocity', poles='south')
    pl_south.camera_position = [(0, 0.01, -3.2), (0, 0, -1), (0, 1, 0)]
    pl_south.screenshot(str(Path(__file__).parent / 'pole_asymmetric_paths_3d_south_pole.png'))
    print('wrote pole_asymmetric_paths_3d_south_pole.png')

    pl_interactive = new_plotter('psi = cos(theta)(1+sin(theta))cos(lambda): velocity field and integration paths')
    pl_interactive.camera_position = [(2.6, 2.6, 1.8), (0, 0, 0), (0, 0, 1)]
    pl_interactive.export_html(str(Path(__file__).parent / 'pole_asymmetric_paths_3d.html'))
    print('wrote pole_asymmetric_paths_3d.html')


if __name__ == '__main__':
    main()
