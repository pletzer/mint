# MINT - Mimetic INTerpolation on the Sphere

<p align="center">
<a href="https://cirrus-ci.com/github/pletzer/mint">
<img src="https://api.cirrus-ci.com/github/pletzer/mint.svg?branch=master"
     alt="Cirrus-CI" /></a>
<a href="https://github.com/pletzer/mint/graphs/contributors">
<img src="https://img.shields.io/github/contributors/pletzer/mint.svg"
     alt="# contributors" /></a>
<a href="https://github.com/pletzer/mint/releases">
<img src="https://img.shields.io/github/v/release/pletzer/mint"
     alt="latest release" /></a>
<a href="https://github.com/pletzer/mint/commits/master">
<img src="https://img.shields.io/github/commits-since/pletzer/mint/latest.svg"
     alt="Commits since last release" /></a>
</p>

----

## Overview

This project contains code to regrid edge centred fields from a source to a destination grid. The grid is stored as a collection of 
grid cells, which have four vertices (i.e. the cells are quadrilaterals). The edge field is stored as integrals of a vector field 
along each edge. The vertex coordinates are stored in longitude-latitude space.

The regridding method is mimetic in the sense that Stokes's theorem is satisfied to near machine precision. In particular, the 
loop integrals of an interpolated vector field deriving from a gradient is zero. 

## Prerequisites

To build the `MINT` Python module:

 * C++ compiler with C++11 support
 * Cython
 * NetCDF library
 * NumPy
 * VTK

We recommend installing the above packages using `conda`.

To build the `MINT` C++ library and the tools:

 * C++ compiler with C++11 support
 * Fortran compiler (e.g., gfortran 6.4)
 * NetCDF library
 * Doxygen
 * VTK

## How to Build the MINT Python Module

The `MINT` Python interface requires `VTK`, `netCDF` and the `tbb` libraries to 
be installed. This is most easily done in a conda environment:
```
conda env create --file requirements/mint.yml
conda activate mint-dev
```

In the root `MINT` directory then type:
```
pip install --no-deps --editable .
```

Check that you can import the mint module:
```
python -c "import mint"
```

To run the tests type:
```
pytest
```
 
## How to Build the MINT C++ Library

Perform the following in order to call `MINT` from `Fortran`, `C` or `C++`:
```
mkdir build
cd build
cmake ..
make -j 8
```

You can specify the compiler with:
```
FC=mpif90 CXX=mpicxx cmake ..; make -j 8
```

You can check that the build was successful by typing:
```
make test
```

### Binary Tools

The above `CMake` build will compile a number of tools. To run the tools, set
`MINT_SRC_DIR` to the location of the `MINT` source directory (e.g. `export MINT_SRC_DIR=..`).

 1. Compute the interpolation weights from a lat-lon grid to the cubed sphere:
 ```
 ./tools/regrid_edges -s $MINT_SRC_DIR/data/latlon4x2.nc:latlon -S 0 \
                      -d $MINT_SRC_DIR/data/cs_4.nc:physics -D 1 \
                      -w weights.nc
 ```

 2. Regrid field `edge_integrated_velocity` from lat-lon to cubed-sphere by loading the previously generated weights:
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


## API Documentation

This is useful if you would like to call `MINT` from `C`, `Python` or `Fortran`:

https://pletzer.github.io/mint/html/


## Conservation Error

The plot below shows the error of mimetic regridding error obtained by computing the 
closed line integral for each cell. 
![alt Error of mimetic regridding](https://raw.githubusercontent.com/pletzer/mint/master/figures/regrid_edgesError.png)
