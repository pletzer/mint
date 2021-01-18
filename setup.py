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
    name="mint", # Replace with your own username
    version=VERSION,
    author_email="alexander@gokliya.net",
    description="Mimetic INterpolation on the sphere",
    url="https://github.com/pletzer/mint",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: ISC License",
    ],
    packages=['mint'],
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
    package_dir={'mint': 'pymint'},
    install_requires=['numpy', 'vtk>=8.1.0', 'netcdf4', 'tbb'],
)
