"""Type of the scannls."""

from typing import Any, Protocol

EventType = tuple[
    str,
    int,
    int,
    tuple[str, str, int, int],
    tuple[int, int, Any],
    tuple[int, int, Any],
    tuple[str, str],
    tuple[str, str],
    list[str],
]


class LoggerType(Protocol):
    """Logger type."""

    def trace(self, msg: str) -> None:
        """Trace."""

    def debug(self, msg: str) -> None:
        """Debug."""

    def info(self, msg: str) -> None:
        """Info."""

    def warning(self, msg: str) -> None:
        """Warning."""

    def error(self, msg: str) -> None:
        """Error."""

    def critical(self, msr: str) -> None:
        """Critical."""

    def success(self, msg: str) -> None:
        """Success."""

    def complete(self) -> Any:
        """Complete."""
