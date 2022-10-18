import copy

from iris.cube import Cube
from iris.coords import DimCoord
import numpy as np


class _DummyMintRegridder:
    def __init__(self, src_coords, tgt_coords, **kwargs):
        self.shape = tgt_coords[0].shape

    def regrid(self, data, dims, **kwargs):
        new_shape = list(data.shape)
        for dim, size in zip(dims, self.shape):
            new_shape[dim] = size
        return np.zeros(self.shape)


def _get_dims(cube):
    # This function may be made more sophisticated to accomodate
    # data not located on a mesh.
    dims = [cube.mesh_dim()]
    return dims


def _get_coords(cube):
    coords = (cube.coord(axis="X"), cube.coord(axis="Y"))
    return coords


def _make_mint_regridder(src_coords, tgt_coords, **kwargs):
    return _DummyMintRegridder(src_coords, tgt_coords, **kwargs)


def _regrid(data, dims, mint_regridder):

    # TODO: replace dummy operation with actual regridding.
    #  Actual operation should make use of mint_regridder.
    new_data = mint_regridder.regrid(data, dims)

    return new_data


def _create_cube(data, src_cube, src_dims, tgt_coords, num_tgt_dims):
    # This function is copied from iris-esmf-regrid.

    new_cube = Cube(data)

    if len(src_dims) == 2:
        grid_dim_x, grid_dim_y = src_dims
    elif len(src_dims) == 1:
        grid_dim_y = src_dims[0]
        grid_dim_x = grid_dim_y + 1
    else:
        raise ValueError(
            f"Source grid must be described by 1 or 2 dimensions, got {len(src_dims)}"
        )
    if num_tgt_dims == 1:
        grid_dim_x = grid_dim_y = min(src_dims)
    for tgt_coord, dim in zip(tgt_coords, (grid_dim_x, grid_dim_y)):
        if len(tgt_coord.shape) == 1:
            if isinstance(tgt_coord, DimCoord):
                new_cube.add_dim_coord(tgt_coord, dim)
            else:
                new_cube.add_aux_coord(tgt_coord, dim)
        else:
            new_cube.add_aux_coord(tgt_coord, (grid_dim_y, grid_dim_x))

    new_cube.metadata = copy.deepcopy(src_cube.metadata)

    def copy_coords(src_coords, add_method):
        for coord in src_coords:
            dims = src_cube.coord_dims(coord)
            if set(src_dims).intersection(set(dims)):
                continue
            offset = num_tgt_dims - len(src_dims)
            dims = [dim if dim < max(src_dims) else dim + offset for dim in dims]
            result_coord = coord.copy()
            # Add result_coord to the owner of add_method.
            add_method(result_coord, dims)

    copy_coords(src_cube.dim_coords, new_cube.add_dim_coord)
    copy_coords(src_cube.aux_coords, new_cube.add_aux_coord)

    return new_cube


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
    # TODO: properly handle the dimensionality of the target
    num_tgt_dims = 1
    return _create_cube(data, cube, dims, tgt_coords, num_tgt_dims)


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
