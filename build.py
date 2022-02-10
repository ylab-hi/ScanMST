# !/usr/bin/env python
"""Build cpp extension.

@Filename:    build.py
@license:     MIT Licence
@Time:        1/7/22 3:00 PM
"""
from pybind11.setup_helpers import build_ext
from pybind11.setup_helpers import Pybind11Extension


def build(setup_kwargs):
    """Build cpp extension."""
    ext_modules = [
        Pybind11Extension(
            "scannls.cppext",
            sources=["src/scannls/cppext/bam.cpp", "src/scannls/cppext/binding.cpp"],
            language="c++",
        )
    ]
    setup_kwargs.update(
        {
            "ext_modules": ext_modules,
            "cmd_class": {"build_ext": build_ext},
            "zip_safe": False,
        }
    )
