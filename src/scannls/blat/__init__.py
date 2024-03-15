"""Blat executable module.

@Filename:    __init__.py
@Author:      YangyangLi
@Time:        5/23/22 10:35 AM
@source:      https://hgdownload.soe.ucsc.edu/admin/exe/
"""

import platform
import sys
from pathlib import Path

from scannls import __PACKAGE_NAME__


def load_blat() -> Path:
    """Load blat.

    @return: Path object.
    """
    blat_path = Path(sys.modules[__PACKAGE_NAME__].__file__).parent / "blat"

    system = platform.system()
    if system == "Windows":
        msg = "Windows is not supported for blat."
        raise NotImplementedError(msg)
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

    Returns
    -------
        fa2bit: fa2bit object.
    """
    path = load_blat() / "faToTwoBit"
    path.chmod(0o755)
    return path
