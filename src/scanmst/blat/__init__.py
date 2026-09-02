"""Blat executable module."""

import os
import platform
import stat
import sys
from pathlib import Path

import requests

from scanmst import __PACKAGE_NAME__
from scanmst.exception import ToolNotFoundError


def get_platform_base_url() -> tuple[str, str]:
    """
    Determine system platform and return the directory name and the
    BASE download URL (root of the architecture folder).

    Returns:
        (system_name, base_url)
    """
    system = platform.system()
    machine = platform.machine()

    if system == "Linux":
        # Base URL for Linux x86_64 binaries. Using an older release because the latest BLAT requires GLIBC 2.33 or higher.
        return "linux", "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/"

    if system == "Darwin":
        if machine == "arm64":
            # Base URL for macOS Apple Silicon
            return "darwin", "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/"
        # Base URL for macOS Intel
        return "darwin", "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/"

    msg = f"Operating system {system} ({machine}) is not supported for blat."
    raise NotImplementedError(msg)


def load_blat() -> Path:
    """Load blat directory path.

    @return: Path object.
    """
    package_file = sys.modules[__PACKAGE_NAME__].__file__
    if package_file is None:
        msg = f"Cannot locate {__PACKAGE_NAME__} on disk"
        raise ToolNotFoundError(msg)
    blat_path = Path(package_file).parent / "blat"
    system_name, _ = get_platform_base_url()
    return blat_path / system_name


def download_blat_tools(logger) -> None:
    """
    Check if blat tools exist for the current platform.
    If not, download them from UCSC.
    """
    try:
        system_name, base_url = get_platform_base_url()
    except NotImplementedError as e:
        logger.error(str(e))
        sys.exit(1)

    # Target directory: src/scannls/blat/{linux|darwin}/
    target_dir = Path(__file__).parent / system_name
    target_dir.mkdir(parents=True, exist_ok=True)

    # Map tool filenames to their relative path on the UCSC server
    # 'blat', 'gfServer', 'gfClient' are inside the 'blat/' subdirectory
    # 'faToTwoBit' is in the root of the architecture directory
    tools_map = {"blat": "blat/blat", "gfServer": "blat/gfServer", "gfClient": "blat/gfClient", "faToTwoBit": "faToTwoBit"}

    for tool_name, relative_path in tools_map.items():
        tool_path = target_dir / tool_name

        # 1. Check if tool exists
        if tool_path.exists():
            # Ensure it is executable
            st = os.stat(tool_path)
            Path(tool_path).chmod(st.st_mode | stat.S_IEXEC)
            continue

        # 2. Construct specific URL
        url = f"{base_url}{relative_path}"
        logger.info(f"Downloading {tool_name} from {url}...")

        try:
            with requests.get(url, stream=True, timeout=60) as r:
                r.raise_for_status()
                with open(tool_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

            # 3. Make executable
            st = os.stat(tool_path)
            Path(tool_path).chmod(st.st_mode | stat.S_IEXEC)
            logger.info(f"Successfully downloaded: {tool_name}")

        except Exception as e:
            logger.error(f"Failed to download {tool_name}: {e}")
            # Clean up partial file
            if tool_path.exists():
                tool_path.unlink()
            sys.exit(f"Error: Could not download {tool_name}. Please check internet connection.")


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
