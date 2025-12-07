"""Main function for scanmst."""

from .cli.arg import parse_args
from .cli.cli import cli


def main():
    """Main function for scanmst."""
    parser = parse_args()
    options = parser.parse_args()

    cli(options)


if __name__ == "__main__":
    main()
