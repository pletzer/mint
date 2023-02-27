import numpy as np
from pathlib import Path
import mint
from iris_utils import _u_v_cubes_from_ugrid_file, _set_vector_field_from_streamfct

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

def test_cs_zt():    

    src_u, src_v = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs128_wind_zt.nc')

    # w2
    _set_vector_field_from_streamfct(src_u, src_v)

    # grid options for a cubed-sphere grid
    src_flags = (1, 1, 1)

    xy = [(0., 0.), (360., 0.), (160., 70.), (0., 0.)]
    flx_calc = mint.IrisMintFlux(src_u.mesh, src_flags=src_flags, tgt_line=xy)
    flux = flx_calc.evaluate_from_vector(src_u, src_v, fs=mint.FUNC_SPACE_W2)

    # check
    assert np.all(np.fabs(flux - 0.0) < 1.e-10)


def test_cs_zt2():    

    src_u, src_v = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs128_wind_zt.nc')

    # w2
    _set_vector_field_from_streamfct(src_u, src_v)

    # grid options for a cubed-sphere grid
    src_flags = (1, 1, 1)

    xybeg = np.array([0., 0.])
    xyend = np.array([256., 72.])
    xy = [xybeg, xyend]
    flx_calc = mint.IrisMintFlux(src_u.mesh, src_flags=src_flags, tgt_line=xy)
    flux = flx_calc.evaluate_from_vector(src_u, src_v, fs=mint.FUNC_SPACE_W2)

    # check
    sbeg = np.cos(xybeg[0]*np.pi/180.)*np.cos(xybeg[1]*np.pi/180.)
    send = np.cos(xyend[0]*np.pi/180.)*np.cos(xyend[1]*np.pi/180.)
    exact = send - sbeg
    assert np.all(np.fabs(flux - exact) < 6.e-5)



