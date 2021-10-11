"""Command-line interface."""
import click


@click.command()
@click.version_option()
def main() -> None:
    """ScanNLS."""


if __name__ == "__main__":
    main(prog_name="ScanNLS")  # pragma: no cover
