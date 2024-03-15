"""Sphinx configuration."""

from datetime import datetime

project = "ScanNLS"
author = "Yangyang Li, Ting-You Wang"
copyright_ = f"{datetime.now().year}, {author}"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinxarg.ext",
    "myst_parser",
]
autodoc_typehints = "description"
html_theme = "furo"
