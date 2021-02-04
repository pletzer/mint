# -*- python -*-
import setuptools
from pathlib import Path
import os
import glob
import re

# extract the version number from the version.txt file
with open("version.txt") as f:
    VERSION = f.read().strip()


# generate pymint/__init__.py from pymint/__init__.py.in
init_file = ""
with open("pymint/__init__.py.in") as fi:
    init_file = re.sub(r'@VERSION@', VERSION, fi.read())
    with open("pymint/__init__.py", 'w') as fo:
        fo.write(init_file)

#
# check the dependencies
#

# VTK
try:
    import vtk
    VTK_VERSION = f'{vtk.VTK_MAJOR_VERSION}.{vtk.VTK_MINOR_VERSION}'
    VTK_INCLUDE_DIR = Path(vtk.__path__[0] + f'/../../../../include/vtk-{VTK_VERSION}')
    VTK_LIBRARIES_DIR = Path(vtk.__path__[0] + '/../../../../lib')
    VTK_LIBRARIES = [f'{lib}-{VTK_VERSION}' for lib in ('vtkCommonComputationalGeometry',
                                                        'vtkIOCore', 
                                                        'vtkIOLegacy', 
                                                        'vtkCommonExecutionModel', 
                                                        'vtkCommonDataModel', 
                                                        'vtkCommonTransforms',
                                                        'vtkCommonMisc',
                                                        'vtkCommonMath',
                                                        'vtkCommonSystem',
                                                        'vtkCommonCore',
                                                        'vtksys')]
except:
    raise RuntimeError('ERROR: "import vtk", must have vtk installed!')

# netCDF
try:
    import netCDF4
    NETCDF_INCLUDE_DIR = Path(netCDF4.__path__[0] + '/../../../../include')
    NETCDF_LIBRARIES_DIR = Path(netCDF4.__path__[0] + '/../../../../lib')
    NETCDF_LIBRARIES = ['netcdf', 'hdf5']
except:
    raise RuntimeError('ERROR: "import netCDF4", must have netcdf4 installed!')

# C++ 11 flag
cpp11_flag = '-std=c++11'
# give a chance to override the C++ 11 flag
cpp_flags = os.getenv("CPPFLAGS")
if cpp_flags:
    cpp11_flag = cpp_flags # on Windows: '/std:c11'


print(f'VTK_VERSION          = {VTK_VERSION}')
print(f'VTK_INCLUDE_DIR      = {VTK_INCLUDE_DIR}')
print(f'VTK_LIBRARIES_DIR    = {VTK_LIBRARIES_DIR}')
print(f'VTK_LIBRARIES        = {VTK_LIBRARIES}')
print(f'NETCDF_INCLUDE_DIR   = {NETCDF_INCLUDE_DIR}')
print(f'NETCDF_LIBRARIES_DIR = {NETCDF_LIBRARIES_DIR}')
print(f'NETCDF_LIBRARIES     = {NETCDF_LIBRARIES}')
print(f'C++11 flag           : {cpp11_flag}')


setuptools.setup(
    name="mintregrid",
    version=VERSION,
    author="Alexander Pletzer",
    author_email="alexander.pletzer@nesi.org.nz",
    description="Mimetic INterpolation on the sphere",
    long_description="""
Interpolates or regrids an edge centred field from a source grid to either a destination grid or a target element. The 
interpolation is mimetic in the sense that line integrals are conserved from source to destination grids, i.e. Stokes'
theorem is statisfied to near machine precision.
    """,
    long_description_content_type="text/x-rst",
    url="https://github.com/pletzer/mint",
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    packages=['mintregrid'],
    ext_modules = [setuptools.Extension('mint', # name of the shared library
                   sources=glob.glob('src/*.cpp'),
                   define_macros=[],
                   include_dirs=['src/',
                                 VTK_INCLUDE_DIR,
                                 NETCDF_INCLUDE_DIR],
                   libraries=VTK_LIBRARIES + NETCDF_LIBRARIES,
                   library_dirs=[vtk.__path__[0]],
                   extra_compile_args=[cpp11_flag,],
                   ),],
    include_package_data=True,
    package_dir={'mintregrid': 'pymint'},
    install_requires=['numpy', 'vtk>=8.1.0', 'netcdf4', 'tbb'],
)
