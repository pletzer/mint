# -*- coding: utf-8 -*-

import glob
import os
from pathlib import Path
import re
import sysconfig
from setuptools import setup, Extension
import sys

from Cython.Build import cythonize

# the following allows us to build a shared library on Windows named libmint.so, see
# https://github.com/himbeles/ctypes-example.git. Works if one uses ctypes to call the
# shared library.
# https://stackoverflow.com/questions/4529555/building-a-ctypes-based-c-library-with-distutils
from distutils.command.build_ext import build_ext as build_ext_orig
class CTypesExtension(Extension): pass
class build_ext(build_ext_orig):

    def build_extension(self, ext):
        self._ctypes = isinstance(ext, CTypesExtension)
        return super().build_extension(ext)

    def get_export_symbols(self, ext):
        if self._ctypes:
            return ext.export_symbols
        return super().get_export_symbols(ext)

    def get_ext_filename(self, ext_name):
        if self._ctypes:
            return ext_name + '.so'
        return super().get_ext_filename(ext_name)


PACKAGE = "mint"


VTK_LIB_NAMES = [
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


def _multiarch_lib_dirs():
    """Extra library search directories used by Debian/Ubuntu-style
    multiarch system installs (e.g. /usr/lib/x86_64-linux-gnu)."""
    dirs = []
    multiarch = sysconfig.get_config_var("MULTIARCH")
    if multiarch:
        dirs.append(Path("/usr/lib") / multiarch)
    dirs += [Path("/usr/lib64"), Path("/usr/lib"), Path("/usr/local/lib"), Path("/usr/local/lib64")]
    return [d for d in dirs if d.is_dir()]


def _lib_dir_has_vtk(lib_dir, version):
    if lib_dir is None or not Path(lib_dir).is_dir():
        return False
    return len(list(Path(lib_dir).glob(f"*vtkIOCore-{version}.*"))) > 0


def getVTK():
    """Locate the VTK headers and libraries needed to build the mint
    extension.

    Resolution order:
      1. Explicit override via the VTK_INCLUDE_DIR / VTK_LIBRARIES_DIR /
         (optionally) VTK_VERSION environment variables -- use this if VTK
         was built from source or installed somewhere non-standard.
      2. A conda/virtualenv environment (sys.exec_prefix/include/vtk-*).
      3. A system package install, e.g. Debian/Ubuntu's "libvtk*-dev",
         Homebrew, or a "make install" under /usr/local.
    """

    env_include = os.getenv("VTK_INCLUDE_DIR")
    env_libdir = os.getenv("VTK_LIBRARIES_DIR")
    if env_include and env_libdir:
        version = os.getenv("VTK_VERSION")
        if not version:
            m = re.search(r"vtk-([0-9]+\.[0-9]+)", env_include)
            if not m:
                raise RuntimeError(
                    "Could not infer VTK_VERSION from VTK_INCLUDE_DIR "
                    f'("{env_include}"); set the VTK_VERSION environment '
                    "variable explicitly (e.g. VTK_VERSION=9.3)."
                )
            version = m.group(1)
        return {
            "VTK_VERSION": version,
            "VTK_INCLUDE_DIR": env_include,
            "VTK_LIBRARIES_DIR": env_libdir,
            "VTK_LIBRARIES": [f"{lib}-{version}" for lib in VTK_LIB_NAMES],
        }

    # candidate (include_root, preferred_lib_dir) pairs to probe, in order
    if os.name == 'nt':
        candidates = [
            (Path(sys.exec_prefix) / "Library" / "include", Path(sys.exec_prefix) / "Library" / "lib"),
        ]
    else:
        candidates = [
            (Path(sys.exec_prefix) / "include", Path(sys.exec_prefix) / "lib"),
            (Path("/usr/local/include"), Path("/usr/local/lib")),
            (Path("/usr/include"), Path("/usr/lib")),
            (Path("/opt/homebrew/include"), Path("/opt/homebrew/lib")),
        ]

    for include_root, preferred_lib_dir in candidates:
        if not include_root.is_dir():
            continue
        matches = sorted(include_root.glob("vtk-*"))
        if not matches:
            continue
        include_dir = matches[-1]
        version = re.sub(r"vtk-", "", include_dir.name)

        lib_dir_candidates = [preferred_lib_dir] + _multiarch_lib_dirs()
        libraries_dir = next(
            (d for d in lib_dir_candidates if _lib_dir_has_vtk(d, version)),
            None,
        )
        if libraries_dir is None:
            # headers found but no matching libraries in any known location;
            # keep looking in case another candidate root matches fully
            continue

        return {
            "VTK_VERSION": version,
            "VTK_INCLUDE_DIR": str(include_dir),
            "VTK_LIBRARIES_DIR": str(libraries_dir),
            "VTK_LIBRARIES": [f"{lib}-{version}" for lib in VTK_LIB_NAMES],
        }

    raise RuntimeError(
        "ERROR: could not locate a VTK development install.\n"
        'Either "conda install -c conda-forge vtk", install your system\'s '
        'VTK development package (e.g. "apt install libvtk9-dev" on '
        "Debian/Ubuntu, \"brew install vtk\" on macOS), build VTK from "
        "source, or point at an existing install by setting the "
        "VTK_INCLUDE_DIR and VTK_LIBRARIES_DIR environment variables "
        "(and VTK_VERSION if it can't be inferred from the include "
        "directory name)."
    )


def getNetCDF():
    """Locate the NetCDF headers and libraries.

    Resolution order:
      1. Explicit override via the NETCDF_INCLUDE_DIR / NETCDF_LIBRARIES_DIR
         environment variables.
      2. "nc-config", the standard way to discover a NetCDF-C install (used
         by the plain CMake build too) -- this works for conda, system
         packages (apt/dnf/brew), and installs built from source alike.
      3. A conda/virtualenv environment or common system prefixes, as a
         last-ditch fallback if nc-config isn't on PATH.
    """
    import shutil
    import subprocess

    env_include = os.getenv("NETCDF_INCLUDE_DIR")
    env_libdir = os.getenv("NETCDF_LIBRARIES_DIR")
    if env_include and env_libdir:
        return {
            "NETCDF_INCLUDE_DIR": env_include,
            "NETCDF_LIBRARIES_DIR": env_libdir,
            "NETCDF_LIBRARIES": ["netcdf", "hdf5"],
        }

    if shutil.which("nc-config"):
        include_dir = subprocess.check_output(
            ["nc-config", "--includedir"], text=True
        ).strip()
        libs_flags = subprocess.check_output(
            ["nc-config", "--libs"], text=True
        ).strip().split()
        libraries_dir = None
        libraries = []
        for flag in libs_flags:
            if flag.startswith("-L"):
                libraries_dir = flag[2:]
            elif flag.startswith("-l"):
                libraries.append(flag[2:])
        if include_dir and libraries_dir and libraries:
            return {
                "NETCDF_INCLUDE_DIR": include_dir,
                "NETCDF_LIBRARIES_DIR": libraries_dir,
                "NETCDF_LIBRARIES": libraries,
            }

    # fall back to probing common prefixes directly
    if os.name == 'nt':
        candidates = [
            (Path(sys.exec_prefix) / "Library" / "include", Path(sys.exec_prefix) / "Library" / "lib"),
        ]
    else:
        candidates = [
            (Path(sys.exec_prefix) / "include", Path(sys.exec_prefix) / "lib"),
            (Path("/usr/local/include"), Path("/usr/local/lib")),
            (Path("/usr/include"), Path("/usr/lib")),
            (Path("/opt/homebrew/include"), Path("/opt/homebrew/lib")),
        ]
    for include_dir, libraries_dir in candidates:
        if (include_dir / "netcdf.h").exists():
            return {
                "NETCDF_INCLUDE_DIR": str(include_dir),
                "NETCDF_LIBRARIES_DIR": str(libraries_dir),
                "NETCDF_LIBRARIES": ["netcdf", "hdf5"],
            }

    raise RuntimeError(
        "ERROR: could not locate a NetCDF-C development install.\n"
        'Either "conda install -c conda-forge libnetcdf", install your '
        'system\'s NetCDF development package (e.g. "apt install '
        'libnetcdf-dev" on Debian/Ubuntu, "brew install netcdf" on macOS), '
        "make sure nc-config is on your PATH, or set the "
        "NETCDF_INCLUDE_DIR and NETCDF_LIBRARIES_DIR environment variables."
    )


# extract the MINT version from file version.txt
with open("version.txt") as f:
    VERSION = f.read().strip()

# generate mint/__init__.py from mint/__init__.py.in
init_file = ""
with open(f"{PACKAGE}/__init__.py.in") as fi:
    init_file = re.sub(r"@VERSION@", VERSION, fi.read())
    with open(f"{PACKAGE}/__init__.py", "w") as fo:
        fo.write(init_file)

vtklib = getVTK()
nclib = getNetCDF()

extra_compile_args = []
cpp_flags = os.getenv("CPPFLAGS")
cxx_flags = os.getenv("CXXFLAGS")
if cxx_flags is not None:
    extra_compile_args = cxx_flags.split()
elif cpp_flags is not None:
    extra_compile_args = cpp_flags.split()

if os.name != 'nt':
    extra_compile_args.append("-std=c++17")

print(f'VTK_VERSION          = {vtklib["VTK_VERSION"]}')
print(f'VTK_INCLUDE_DIR      = {vtklib["VTK_INCLUDE_DIR"]}')
print(f'VTK_LIBRARIES_DIR    = {vtklib["VTK_LIBRARIES_DIR"]}')
print(f'VTK_LIBRARIES        = {vtklib["VTK_LIBRARIES"]}')
print(f'NETCDF_INCLUDE_DIR   = {nclib["NETCDF_INCLUDE_DIR"]}')
print(f'NETCDF_LIBRARIES_DIR = {nclib["NETCDF_LIBRARIES_DIR"]}')
print(f'NETCDF_LIBRARIES     = {nclib["NETCDF_LIBRARIES"]}')
print(f"extra_compile_args   = {extra_compile_args}")

extensions = [
    CTypesExtension(
        f"lib{PACKAGE}",
        sources=glob.glob("src/*.cpp"),
        define_macros=[],
        include_dirs=["src/", vtklib["VTK_INCLUDE_DIR"], nclib["NETCDF_INCLUDE_DIR"]],
        libraries=vtklib["VTK_LIBRARIES"] + nclib["NETCDF_LIBRARIES"],
        library_dirs=[vtklib["VTK_LIBRARIES_DIR"], nclib["NETCDF_LIBRARIES_DIR"]],
        extra_compile_args=extra_compile_args,
        language="c++",
    )
]

setup(
    ext_modules=cythonize(extensions, compiler_directives=dict(language_level=3)),
    cmdclass={'build_ext': build_ext},
)
