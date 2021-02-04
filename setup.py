# -*- python -*-
import setuptools
from pathlib import Path
import os
import glob
import re
import sys


def getVtk():
    inc_dir = Path(sys.exec_prefix) / Path('include')
    vtk_version = list(inc_dir.glob('vtk-*'))[-1].name
    vtk_version = re.sub(r'vtk-', '', vtk_version)
    return {'VTK_VERSION': vtk_version,
            'VTK_INCLUDE_DIR': str( Path(sys.exec_prefix) / Path('include') / Path(f'vtk-{vtk_version}') ),
            'VTK_LIBRARIES_DIR': str( Path(sys.exec_prefix) / Path('lib') )}


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
# VTK installation. Location  can be set by the user via environment variables, if
# desired. Otherwise will infer form a typical pip/conda installation
#
vtk_libraries = ('vtkCommonComputationalGeometry', 'vtkIOCore', 
                 'vtkIOLegacy', 'vtkCommonExecutionModel', 'vtkCommonDataModel', 
                 'vtkCommonTransforms', 'vtkCommonMisc', 'vtkCommonMath', 'vtkCommonSystem',
                 'vtkCommonCore', 'vtksys')

VTK_VERSION = os.getenv('VTK_VERSION')
VTK_INCLUDE_DIR = os.getenv('VTK_INCLUDE_DIR')
VTK_LIBRARIES_DIR = os.getenv('VTK_LIBRARIES_DIR')
if not (VTK_VERSION and VTK_INCLUDE_DIR and VTK_LIBRARIES_DIR):
    lib = getVtk()
    VTK_VERSION = lib['VTK_VERSION']
    VTK_INCLUDE_DIR = lib['VTK_INCLUDE_DIR']
    VTK_LIBRARIES_DIR = lib['VTK_LIBRARIES_DIR']
VTK_LIBRARIES = [f'{lib}-{VTK_VERSION}' for lib in vtk_libraries]

#
# netCDF installation. 
#
NETCDF_INCLUDE_DIR = os.getenv('NETCDF_INCLUDE_DIR')
NETCDF_LIBRARIES_DIR = os.getenv('NETCDF_LIBRARIES_DIR')
if not (NETCDF_INCLUDE_DIR and NETCDF_LIBRARIES_DIR):
    NETCDF_INCLUDE_DIR = str( Path(sys.exec_prefix) / Path('include') )
    NETCDF_LIBRARIES_DIR = str( Path(sys.exec_prefix) / Path('lib') )
NETCDF_LIBRARIES = ['netcdf', 'hdf5']

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
                   library_dirs=[VTK_INCLUDE_DIR, NETCDF_LIBRARIES_DIR],
                   extra_compile_args=[cpp11_flag,],
                   ),],
    include_package_data=True,
    package_dir={'mintregrid': 'pymint'},
    install_requires=['numpy', 'vtk>=8.1.0', 'netcdf4', 'tbb'],
)
