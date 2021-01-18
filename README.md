MINT - Mimetic Interpolation on the Sphere

## Overview

This project contains code to regrid edge centred fields from a source to a destination grid. The grid is stored as a collection of 
grid cells, which have four vertices (i.e. the cells are quadrilaterals). The edge field is stored as integrals of a vector field 
along each edge. The vertex coordinates are stored in longitude-latitude space.

The regridding method is mimetic in the sense that Stokes's theorem is satisfied to near machine precision. In particular, the 
loop integrals of an interpolated vector field deriving from a gradient is zero. 

## Prerequisites

You will need for the python module:

 * C++ compiler with C++ 11 support
 * netcdf4
 * vtk
 * numpy

We recommend to install the above packages with Anaconda.

For building the library and the tools

 * C++ compiler with C++ 11 support
 * Fortran compiler (e.g gfortran 6.4)
 * NetCDF library
 * Doxygen


## How to build the python mint module

```
python setup.py install 
```

Check that you can import the mint module
```
python -c "import mint"
```

To run a quick test type
```
python pymint/tests/test_regrid_edges.py
```
 
## How to build the MINT library

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

## API Documentation

This is useful if you would like to call MINT from C, Python or Fortran:

https://pletzer.github.io/mint/html/

## Examples

Set MINT_SRC_DIR to the location of the MINT source directory (e.g. `export MINT_SRC_DIR=..`).

 1. Compute the interpolation weights from a lat-lon grid to the cubed sphere:
 ```
 ./tools/regrid_edges -s $MINT_SRC_DIR/data/latlon4x2.nc:latlon -S 0 \
                      -d $MINT_SRC_DIR/data/cs_4.nc:physics -D 1 \
                      -w weights.nc
 ```

 2. Regrid field "edge_integrated_velocity" from lat-lon to cubed-sphere by loading the previously generated weights:
 ```
 ./tools/regrid_edges -s $MINT_SRC_DIR/data/latlon4x2.nc:latlon -S 0 \
                      -d $MINT_SRC_DIR/data/cs_4.nc:physics -D 1 \
                      -v edge_integrated_velocity -W weights.nc -o edge_integrated_velocity.vtk
 ```

 3. Compute the weights and regrid in one step:
 ```
 ./tools/regrid_edges -s $MINT_SRC_DIR/data/latlon4x2.nc:latlon -S 0 \
                      -d $MINT_SRC_DIR/data/cs_4.nc:physics -D 1 \
                      -v edge_integrated_velocity -W -o edge_integrated_velocity.vtk
 ```

 4. Regrid a time dependent field with elevation:
 ```
./tools/regrid_edges -s $MINT_SRC_DIR/data/lonlatzt_8x4x3x2.nc:mesh2d -S 0 \
                     -d $MINT_SRC_DIR//data/c24_u_integrated.nc:physics -D 1 \
                     -v u@$MINT_SRC_DIR/data/lonlatzt_8x4x3x2.nc \
                     -O regridedges_xyzt.nc -o regridedges_xyzt.vtk

 ```

## Conservation error

Error of mimetic regridding:
![alt Error of mimetic regridding](https://raw.githubusercontent.com/pletzer/mint/master/figures/regrid_edgesError.png)









