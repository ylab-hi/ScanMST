#!/usr/bin/env python3
"""Assert a built wheel is shaped correctly before it can ever be published.

Two failure modes this guards against, both of which produce a wheel that
installs cleanly and then breaks at import time:

* the extension landing somewhere other than ``scanmst/_cppext*.so`` --
  ``scanmst/cppext/__init__.py`` does ``from scanmst._cppext.cppext import *``
* ``libhts`` not being bundled, so the wheel depends on a system htslib the
  user almost certainly does not have

Usage:
    python scripts/check_wheel.py wheelhouse/*.whl
    python scripts/check_wheel.py --no-bundled-libs dist/*.whl
"""

from __future__ import annotations

import argparse
import re
import sys
import zipfile
from pathlib import Path

EXTENSION_RE = re.compile(r"^scanmst/_cppext\.[^/]+\.so$")
JUNK_PATTERNS = ("blat/linux/", "blat/darwin/", "__pycache__/", ".omc/")


def check(path: Path, *, expect_bundled_libs: bool) -> list[str]:
    """Return a list of problems found in the wheel at *path*."""
    problems: list[str] = []
    names = zipfile.ZipFile(path).namelist()

    if not any(EXTENSION_RE.match(n) for n in names):
        misplaced = [n for n in names if "_cppext" in n and n.endswith(".so")]
        problems.append(
            f"extension not at scanmst/_cppext*.so (found: {misplaced or 'nothing'})",
        )

    junk = [n for n in names if any(p in n for p in JUNK_PATTERNS)]
    if junk:
        problems.append(f"junk in wheel: {junk[:5]}")

    bundled = [n for n in names if ".libs/" in n]
    if expect_bundled_libs:
        if not any("libhts" in n for n in bundled):
            problems.append(f"libhts not bundled (auditwheel output: {bundled or 'none'})")
        if "linux_x86_64.whl" in path.name:
            problems.append("wheel is not manylinux-tagged; PyPI will reject it")

    print(f"== {path.name}")
    for n in sorted(names):
        if EXTENSION_RE.match(n) or ".libs/" in n:
            print(f"   {n}")

    return problems


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("wheels", nargs="+", type=Path)
    parser.add_argument(
        "--no-bundled-libs",
        action="store_true",
        help="skip auditwheel-specific checks (for plain local builds)",
    )
    args = parser.parse_args()

    failed = False
    for wheel in args.wheels:
        problems = check(wheel, expect_bundled_libs=not args.no_bundled_libs)
        for problem in problems:
            print(f"   FAIL: {problem}")
            failed = True
        if not problems:
            print("   OK")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
