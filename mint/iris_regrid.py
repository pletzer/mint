from iris.cube import Cube
import numpy as np


def _get_dims(cube):
    # This function may be made more sophisticated to accomodate
    # data not located on a mesh.
    dims = [cube.mesh_dim()]
    return dims


def _get_coords(cube):
    coords = (cube.coord(axis="X"), cube.coord(axis="Y"))
    return coords


def _make_mint_regridder(src_coords, tgt_coords, **kwargs):
    return None


def _regrid(data, dims, regrid_info):
    tgt_coords, mint_regridder = regrid_info

    # TODO: replace dummy operation with actual regridding.
    #  Actual operation should make use of mint_regridder.
    out_shape = list[data.shape]
    out_shape[dims[0]] = tgt_coords[0].shape[0]
    new_data = np.oneslike(out_shape)

    return new_data


def _create_cube(data, src, src_dims, tgt_coords):
    # TODO: Add metadata to new cube.
    result_cube = Cube(data)
    return result_cube


def _prepare(src, tgt, **kwargs):
    src_coords = _get_coords(src)
    tgt_coords = _get_coords(tgt)
    mint_regridder = _make_mint_regridder(src_coords, tgt_coords, **kwargs)
    regrid_info = [tgt_coords, mint_regridder]
    return regrid_info


def _perform(cube, regrid_info):
    tgt_coords, mint_regridder = regrid_info
    dims = _get_dims(cube)
    data = _regrid(cube.data, dims, mint_regridder)
    return _create_cube(data, cube, dims, tgt_coords)


class _MINTRegridder:
    def __init__(self, src, tgt, **kwargs):
        self.regrid_info = _prepare(src, tgt, **kwargs)

    def __call__(self, src):
        return _perform(src, self.regrid_info)


class MINTScheme:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def regridder(self, src, tgt):
        return _MINTRegridder(src, tgt, **self.kwargs)
