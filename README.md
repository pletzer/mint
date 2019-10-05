MINT - Mimetic Interpolation on the Sphere

## Overview

This project contains code to regrid edge centred fields from a source to a destination grid. The grid is stored as a collection of 
grid cells, which have four vertices (i.e. the cells are quadrilaterals). The edge field is stored as integrals of a vector field 
along each edge. The vertex coordinates are stored in longitude-latitude space.

The regridding method is mimetic in the sense that Stokes's theorem is satisfied to near machine precision. In particular, the 
loop integrals of an interpolated vector field deriving from a gradient is zero. 

## Prerequisites

You will need to have:

 * Python (tested with 2.7.14 and 3.6.6)
 * numpy (tested with 1.14.2 and 1.15.4)
 * netCDF4 (tested with 1.3.1 and 1.4.0)
 * VTK with python bindings enabled (tested with 8.1.0 and 8.1.1)


 * C++ compiler (e.g. g++ 6.4)
 * Fortran compiler (e.g gfortran 6.4)
 * NetCDF library


 We recommend to install Python and related packages using Anaconda.
 
## How to build MINT

```
mkdir build
cd build
cmake ..
make -j 8
```

You can specify the compiler with
```
FC=mpif90 CXX=mpicxx cmake ..; make -j 8
```

## Checking that the build was sucessful

```
make test
```

## Example how to regrid from lat-lon to a cubed-sphere grid

In directory `tools/`:
```
cd ../tools
```

 1. Generate lat-lon grid of Unified Model NetCDF flavour
 ```
 python latlon.py -nlon 100 -nlat 60 -o um100x60.nc
 ```

 2. Convert the destination grid to VTK file format
 ```
 python ugrid_reader.py -i ../data/cs_16.nc:physics -V cs_16.vtk
 ```

 3. Read the source grid, generate edge data and save the result as a VTK file
 ```
 python latlon_reader.py -i um100x60.nc -stream "sin(2*pi*x/360.)*cos(pi*y/180.)" -V um100x60.vtk
 ```
 Note: this sets the stream function where x, y are the longitude, respectively, latitudes in degrees. 
 The vector field integral on edges is the difference of stream function values between the edge vertices.

The figure below shows a detail of the source (red) and destination (green) grids.
![alt Source (red) and destination (green) grids](https://raw.githubusercontent.com/pletzer/mint/master/figures/srcAndDstGrids.png)


 4. Regrid the above field from the UM source grid to a cubed-sphere and save the result in a VTK file
 ```
 python regrid_edges.py -s um100x60.vtk -v "edge_integrated_velocity" -d cs_16.vtk -o regrid_edges.vtk
 ```
 or, alternatively, the C++ version:
 ```
 ../build/tools/regrid_edges -s um100x60.vtk -v "edge_integrated_velocity" -d cs_16.vtk -o regrid_edges.vtk
 ```


 The above should print very small values, indicating that the loop integrals are nearly zero for each cell, and thus satisfying Stokes' theorem.

It is important to note that treating the the edge field as a nodal field, i.e. applying bilinear regridding 
does not produce zero loop integrals:
```
python regrid_verts.py -s um100x60.vtk -v "edge_integrated_velocity" -d cs_16.vtk -o regrid_verts.vtk
```

Error of bilinear regridding:
![alt Error of bilinear regridding](https://raw.githubusercontent.com/pletzer/mint/master/figures/regrid_vertsError.png)

Error of mimetic regridding:
![alt Error of mimetic regridding](https://raw.githubusercontent.com/pletzer/mint/master/figures/regrid_edgesError.png)









