import copy

from iris.cube import Cube
from iris.coords import DimCoord, AuxCoord
from iris.experimental.ugrid.mesh import Mesh
import numpy as np

import mint


class IrisToMintMeshAdaptor:

    def __init__(self, iris_mesh: Mesh, flags: tuple):
        """
        Create a MINT grid from an Iris mesh
        :param iris_mesh: Iris mesh object
        :param flags: flags to pass to the grid. Example (0, 0, 1) for a regular grid 
                      and (1, 1, 1) for a cubed-sphere grid. 
        """

        # vertex points
        x = iris_mesh.node_coords.node_x.points
        y = iris_mesh.node_coords.node_y.points
        num_points = x.shape[0]

        # needs to be 3d
        self.points = np.zeros((num_points, 3), np.float64)
        self.points[:, 0] = x
        self.points[:, 1] = y

        # zero-based connecticvity
        self.face2nodes = np.array(iris_mesh.face_node_connectivity.indices_by_location() - \
                                   iris_mesh.face_node_connectivity.start_index, np.uint64)
        self.edge2nodes = np.array(iris_mesh.edge_node_connectivity.indices_by_location() - \
                                   iris_mesh.edge_node_connectivity.start_index, np.uint64)
        
        self.grid = mint.Grid()
        self.grid.setFlags(flags[0], flags[1], flags[2])
        self.grid.loadFromUgrid2DData(self.points, self.face2nodes, self.edge2nodes)

        self.num_faces = self.face2nodes.shape[0]
        self.num_edges = self.edge2nodes.shape[0]
        self.num_points = num_points

    def get_grid(self):
        return self.grid


def createIrisCube(data: np.ndarray, src_cube: src_cube, src_dims: tuple, tgt_coords, num_tgt_dims):
    """
    Create an Iris cube.
    """

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


