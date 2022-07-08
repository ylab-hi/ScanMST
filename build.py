# !/usr/bin/env python
"""Build cpp extension.

@Filename:    build.py
@license:     MIT Licence
@Time:        1/7/22 3:00 PM
"""
import os
from contextlib import contextmanager
from functools import wraps

from pybind11.setup_helpers import build_ext
from pybind11.setup_helpers import Pybind11Extension

HTSLIB_LIBRARY_DIR = os.environ.get("HTSLIB_LIBRARY_DIR", None)
HTSLIB_INCLUDE_DIR = os.environ.get("HTSLIB_INCLUDE_DIR", None)
HTSLIB_CONFIGURE_OPTIONS = os.environ.get("HTSLIB_CONFIGURE_OPTIONS", None)
HTSLIB_SOURCE = None

if not HTSLIB_LIBRARY_DIR:
    raise ValueError("HTSLIB_LIBRARY_DIR is not set")

# linking against a shared, externally installed htslib version, no
# sources required for htslib
htslib_sources = []
shared_htslib_sources = []
chtslib_sources = []
htslib_library_dirs = [HTSLIB_LIBRARY_DIR]
htslib_include_dirs = [HTSLIB_INCLUDE_DIR]
external_htslib_libraries = ["z", "hts"]


@contextmanager
def change_dir(path: str):
    """Change directory."""
    save_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(save_dir)


def change_env(key: str, value: str):
    """Change environment variable."""

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            old_env = os.environ.get(key, None)
            os.environ[key] = old_env + " " + value if old_env else value
            func(*args, **kwargs)
            os.environ[key] = old_env if old_env else " "

        return wrapper

    return decorator


@change_env("CPPFLAGS", "-o2")
def build(setup_kwargs):
    """Build cpp extension."""
    ext_modules = [
        Pybind11Extension(
            "scannls._cppext",
            sources=[
                "src/scannls/cppext/src/bam.cpp",
                "src/scannls/cppext/src/rescuer.cpp",
                "src/scannls/cppext/src/ssw.c",
                "src/scannls/cppext/src/ssw_cpp.cpp",
                "src/scannls/cppext/src/binding.cpp",
            ],
            include_dirs=htslib_include_dirs + ["src/scannls/cppext/include"],
            library_dirs=htslib_library_dirs,
            libraries=external_htslib_libraries,
        )
    ]
    setup_kwargs.update(
        {
            "ext_modules": ext_modules,
            "cmdclass": {"build_ext": build_ext},
            "zip_safe": False,
        }
    )
