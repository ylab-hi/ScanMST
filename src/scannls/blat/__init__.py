"""Blat executable module.

@Filename:    __init__.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        5/23/22 10:35 AM
@source: https://hgdownload.soe.ucsc.edu/admin/exe/
"""
import platform
from importlib import resources
from pathlib import Path

PACKAGE_NAME = "scannls_ont"


def load_blat() -> Path:
    """Load blat.

    @return: Path object.
    """
    with resources.path(PACKAGE_NAME, "blat") as f:
        blat_path = f
    system = platform.system()
    if system == "Windows":
        raise NotImplementedError("Windows is not supported for blat.")
    return blat_path / system.lower()


def load_gfserver() -> Path:
    """Load gfServer.

    @return: Path object.
    """
    path = load_blat() / "gfServer"
    path.chmod(0o755)
    return path


def load_gfclient() -> Path:
    """Load gfClient.

    @return: Path object.
    """
    path = load_blat() / "gfClient"
    path.chmod(0o755)
    return path


def load_fa2bit():
    """Load fa2bit.

    Returns:
        fa2bit: fa2bit object.
    """
    path = load_blat() / "faToTwoBit"
    path.chmod(0o755)
    return path
