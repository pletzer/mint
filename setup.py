# -*- coding: utf-8 -*-

import glob
import os
from pathlib import Path
import re
from setuptools import setup, Extension
import sys

from Cython.Build import cythonize


PACKAGE = "mint"


def getCondaVTK():
    """get the VTK installed by conda"""

    vtk_libs = [
        "vtkCommonComputationalGeometry",
        "vtkIOCore",
        "vtkIOLegacy",
        "vtkCommonExecutionModel",
        "vtkCommonDataModel",
        "vtkCommonTransforms",
        "vtkCommonMisc",
        "vtkCommonMath",
        "vtkCommonSystem",
        "vtkCommonCore",
        "vtksys",
    ]

    include_dir = Path(sys.exec_prefix) / Path("include")

    try:
        version = list(include_dir.glob("vtk-*"))[-1].name
    except IndexError:
        raise RuntimeError('ERROR: you need to "conda install vtk"')

    version = re.sub(r"vtk-", "", version)
    include_dir = str(include_dir / Path(f"vtk-{version}"))
    libraries_dir = str(Path(sys.exec_prefix) / Path("lib"))
    libraries = [f"{lib}-{version}" for lib in vtk_libs]
    result = {
        "VTK_VERSION": version,
        "VTK_INCLUDE_DIR": include_dir,
        "VTK_LIBRARIES_DIR": libraries_dir,
        "VTK_LIBRARIES": libraries,
    }
    return result


def getCondaNetCDF():
    """Get the NetCDF installed by conda."""

    include_dir = str(Path(sys.exec_prefix) / Path("include"))
    libraries_dir = str(Path(sys.exec_prefix) / Path("lib"))
    libraries = ["netcdf", "hdf5"]

    if not (include_dir / Path("netcdf.h")).exists():
        raise RuntimeError('ERROR: you need to "conda install libnetcdf"')

    result = {
        "NETCDF_INCLUDE_DIR": include_dir,
        "NETCDF_LIBRARIES_DIR": libraries_dir,
        "NETCDF_LIBRARIES": libraries,
    }
    return result


# extract the MINT version from file version.txt
with open("version.txt") as f:
    VERSION = f.read().strip()

# generate mint/__init__.py from mint/__init__.py.in
init_file = ""
with open(f"{PACKAGE}/__init__.py.in") as fi:
    init_file = re.sub(r"@VERSION@", VERSION, fi.read())
    with open(f"{PACKAGE}/__init__.py", "w") as fo:
        fo.write(init_file)

vtklib = getCondaVTK()
nclib = getCondaNetCDF()

cpp_flags = os.getenv("CPPFLAGS")
cxx_flags = os.getenv("CXXFLAGS")
if cxx_flags is not None:
    extra_compile_args = cxx_flags.split()
elif cpp_flags is not None:
    extra_compile_args = cpp_flags.split()
else:
    extra_compile_args = ["-std=c++11"]

print(f'VTK_VERSION          = {vtklib["VTK_VERSION"]}')
print(f'VTK_INCLUDE_DIR      = {vtklib["VTK_INCLUDE_DIR"]}')
print(f'VTK_LIBRARIES_DIR    = {vtklib["VTK_LIBRARIES_DIR"]}')
print(f'VTK_LIBRARIES        = {vtklib["VTK_LIBRARIES"]}')
print(f'NETCDF_INCLUDE_DIR   = {nclib["NETCDF_INCLUDE_DIR"]}')
print(f'NETCDF_LIBRARIES_DIR = {nclib["NETCDF_LIBRARIES_DIR"]}')
print(f'NETCDF_LIBRARIES     = {nclib["NETCDF_LIBRARIES"]}')
print(f"extra_compile_args   = {extra_compile_args}")

extensions = [
    Extension(
        f"lib{PACKAGE}",
        sources=glob.glob("src/*.cpp"),
        define_macros=[],
        include_dirs=["src/", vtklib["VTK_INCLUDE_DIR"], nclib["NETCDF_INCLUDE_DIR"]],
        libraries=vtklib["VTK_LIBRARIES"] + nclib["NETCDF_LIBRARIES"],
        library_dirs=[vtklib["VTK_INCLUDE_DIR"], nclib["NETCDF_LIBRARIES_DIR"]],
        extra_compile_args=extra_compile_args,
        language="c++",
    )
]

setup(
    ext_modules=cythonize(extensions, compiler_directives=dict(language_level=3)),
)
