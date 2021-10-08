# MINT - Mimetic INTerpolation on the Sphere

<p align="left">
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

## References

`MINT` implements the 2D+1 interpolation and regridding method described in 
[Mimetic Interpolation of Vector Fields on Arakawa C/D Grids](https://journals.ametsoc.org/view/journals/mwre/147/1/mwr-d-18-0146.1.xml). 
See [Conservative interpolation of edge and face data on n dimensional structured grids using differential forms](https://www.sciencedirect.com/science/article/pii/S0021999115005562?via%3Dihub) for the corresponding methods in n dimensions.


## How to build `MINT` as  a Python module

This the most user friendly access to the MINT package. Not all the functionality from C++ is exposed but it contains the most important parts. (Addition features can be added on request.)

We recommend you install `Miniconda3`. Once you have `Miniconda3` installed, type (Unix)
```
conda env create --file requirements/mint.yml
conda activate mint-dev
```

If you find the above commands to fail on Windows, you can try instead (in the Anaconda prompt terminal):
```
conda create -n mint-dev
conda activate mint-dev
conda install -c conda-forge cmake cython setuptools tbb-devel pip libnetcdf vtk=9.0.3 numpy pytest
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

Documentation of the API can be obtained by typing:
```
pydoc mint
pydoc mint.regrid_edges
pydoc mint.polyline_integral
```
for instance.

 
## How to Build the MINT C++ Library

If you want to call `MINT` from `Fortran`, `C` or `C++`:
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

It is also possible to turn off Fortran code and specify the location of the `VTK` library:
```
cmake -DBUILD_FORTRAN=OFF -DVTK_INCLUDE_DIR=/usr/local/Cellar/vtk/9.0.3/include/vtk-9.0/ -DVTK_LIBRARY_DIR=/usr/local/Cellar/vtk/9.0.3/lib ..
make -j 4
```

On Windows 10 Entreprise, in the Anaconda prompt terminal I do:
```
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
cmake -G "NMake Makefiles" -DBUILD_FORTRAN=OFF -DVTK_INCLUDE_DIR="%userprofile%\miniconda3\envs\mint-dev\Library\include\vtk-9.0" -DVTK_LIBRARY_DIR="%userprofile%\miniconda3\envs\mint-dev\Library\lib" -DNETCDF_INCLUDE_DIRS="%userprofile%\miniconda3\envs\mint-dev\Library\include" -DNETCDF_LIBRARIES="%userprofile%\miniconda3\envs\mint-dev\Library\lib\netcdf.lib" ..
nmake
```

You can check that the build was successful by typing:
```
make test
```


## API Documentation

This is useful if you would like to call `MINT` from `C`, `Python` or `Fortran`:

https://pletzer.github.io/mint/html/


## Examples

### Conservation Error

The plot below shows the error of mimetic regridding error obtained by computing the 
closed line integral for each cell on the cubed-sphere. 
![alt Error of mimetic regridding](https://raw.githubusercontent.com/pletzer/mint/master/figures/regrid_edgesError.png)

### Computing the water flow across irregular boundaries

![alt Gulf Stream](https://raw.githubusercontent.com/pletzer/mint/master/figures/gulfStream.png)
![Around New Zealand](https://raw.githubusercontent.com/pletzer/mint/master/figures/nz.png) 
