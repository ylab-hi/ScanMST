#!/usr/bin/env python3
"""Main function for scannls."""
import sys

from scannls import cli
from scannls.cli import parse_args


def main():
    """Main function for scannls."""
    if sys.version_info < (3, 8):
        sys.exit(
            "Sorry, this code need Python 3.8 or higher. Please update. Aborting..."
        )
    parser = parse_args()

    if len(sys.argv[1:]) < 1:
        parser.print_help()
        raise SystemExit
    else:
        options = parser.parse_args()

    try:
        cli(options)
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)


if __name__ == "__main__":
    main()
