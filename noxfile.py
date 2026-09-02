"""Nox sessions."""

import shutil
import sys
from pathlib import Path
from textwrap import dedent

import nox
from nox import Session

package = "scanmst"
# Capped at 3.10 by the pysam==0.19.0 pin (no cp311/cp312 wheels), matching
# requires-python in pyproject.toml.
python_versions = ["3.9", "3.10"]
nox.needs_version = ">= 2024.3.2"
# uv provides the virtualenvs; fall back to virtualenv where uv is unavailable.
nox.options.default_venv_backend = "uv|virtualenv"
nox.options.sessions = [
    "pre-commit",
    "safety",
    "mypy",
    "tests",
    "typeguard",
]


def activate_virtualenv_in_precommit_hooks(session: Session) -> None:
    """Activate virtualenv in hooks installed by pre-commit.

    This function patches git hooks installed by pre-commit to activate the
    session's virtual environment. This allows pre-commit to locate hooks in
    that environment when invoked from git.

    Args:
    ----
        session: The Session object.
    """
    if session.bin is None:
        return

    virtualenv = session.env.get("VIRTUAL_ENV")
    if virtualenv is None:
        return

    hookdir = Path(".git") / "hooks"
    if not hookdir.is_dir():
        return

    for hook in hookdir.iterdir():
        if hook.name.endswith(".sample") or not hook.is_file():
            continue

        text = hook.read_text()
        bindir = repr(session.bin)[1:-1]  # strip quotes
        if not ((Path("A") == Path("a") and bindir.lower() in text.lower()) or bindir in text):
            continue

        lines = text.splitlines()
        if not (lines[0].startswith("#!") and "python" in lines[0].lower()):
            continue

        header = dedent(
            f"""\
            import os
            os.environ["VIRTUAL_ENV"] = {virtualenv!r}
            os.environ["PATH"] = os.pathsep.join((
                {session.bin!r},
                os.environ.get("PATH", ""),
            ))
            """,
        )

        lines.insert(1, header)
        hook.write_text("\n".join(lines))


@nox.session(name="pre-commit", python=python_versions)
def precommit(session: Session) -> None:
    """Lint using pre-commit."""
    args = session.posargs or ["run", "--all-files", "--show-diff-on-failure"]
    # .pre-commit-config.yaml pins its own hook environments, so only
    # pre-commit itself is needed here.
    session.install("pre-commit", "pre-commit-hooks")
    session.run("pre-commit", *args)
    if args and args[0] == "install":
        activate_virtualenv_in_precommit_hooks(session)


@nox.session(python=python_versions)
def safety(session: Session) -> None:
    """Scan dependencies for known vulnerabilities."""
    requirements = Path(session.create_tmp()) / "requirements.txt"
    # Replaces nox-poetry's session.poetry.export_requirements().
    session.run(
        "uv",
        "export",
        "--format",
        "requirements-txt",
        "--no-emit-project",
        "--no-hashes",
        "-o",
        str(requirements),
        external=True,
    )
    # `safety check` was deprecated in June 2024 and now exits 64 regardless of
    # findings; `safety scan` requires an account and a SAFETY_API_KEY secret.
    # pip-audit needs neither and reads the same requirements file.
    session.install("pip-audit")
    # --no-deps/--disable-pip: `uv export` already emits a fully pinned, resolved
    # set, so there is nothing to re-resolve. Without them pip-audit shells out to
    # pip and tries to build every sdist (pysam among them) just to rediscover
    # dependencies it already has.
    session.run(
        "pip-audit",
        "--no-deps",
        "--disable-pip",
        "--requirement",
        str(requirements),
    )


@nox.session(python=python_versions)
def mypy(session: Session) -> None:
    """Type-check using mypy."""
    # docs/conf.py was removed from this list: the docs are mkdocs, not sphinx,
    # and that file has never existed in this repo.
    args = session.posargs or ["src", "tests"]
    session.install(".")
    session.install("mypy", "pytest")
    session.run("mypy", *args)
    if not session.posargs:
        session.run("mypy", f"--python-executable={sys.executable}", "noxfile.py")


@nox.session(python=python_versions)
def tests(session: Session) -> None:
    """Run the test suite."""
    session.install(".")
    session.install("coverage[toml]", "pytest", "pygments", "pytest-mock")
    try:
        session.run("coverage", "run", "--parallel", "-m", "pytest", *session.posargs)
    finally:
        if session.interactive:
            session.notify("coverage", posargs=[])


@nox.session
def coverage(session: Session) -> None:
    """Produce the coverage report."""
    args = session.posargs or ["report", "-i"]

    session.install("coverage[toml]")

    if not session.posargs and any(Path().glob(".coverage.*")):
        session.run("coverage", "combine")

    session.run("coverage", *args)


@nox.session(python=python_versions)
def typeguard(session: Session) -> None:
    """Runtime type checking using Typeguard."""
    session.install(".")
    session.install("pytest", "typeguard", "pygments", "pytest-mock")
    session.run("pytest", f"--typeguard-packages={package}", *session.posargs)


@nox.session(python="3.10")
def refurb(session: Session) -> None:
    """Refactoring suggestions using refurb."""
    # Previously this passed `args` to session.install instead of session.run,
    # so refurb never actually ran.
    args = session.posargs or ["src", "tests"]
    session.install(".")
    session.install("refurb")
    session.run("refurb", *args)


@nox.session(python=python_versions)
def xdoctest(session: Session) -> None:
    """Run examples with xdoctest."""
    args = session.posargs or ["all"]
    session.install(".")
    session.install("xdoctest[colors]")
    session.run("python", "-m", "xdoctest", package, *args)


@nox.session(name="docs-build", python=python_versions)
def docs_build(session: Session) -> None:
    """Build the documentation."""
    args = session.posargs or ["build", "--strict"]
    session.install(".")
    session.install(
        "mkdocs",
        "mkdocs-material",
        "mkdocs-minify-plugin",
        "mkdocs-git-revision-date-localized-plugin",
        "mkdocstrings[python]",
        "pymdown-extensions",
    )

    build_dir = Path("site")
    if build_dir.exists():
        shutil.rmtree(build_dir)

    session.run("mkdocs", *args)


@nox.session(python=python_versions)
def docs(session: Session) -> None:
    """Build and serve the documentation with live reloading on file changes."""
    args = session.posargs or ["serve"]
    session.install(".")
    session.install(
        "mkdocs",
        "mkdocs-material",
        "mkdocs-minify-plugin",
        "mkdocs-git-revision-date-localized-plugin",
        "mkdocstrings[python]",
        "pymdown-extensions",
    )

    session.run("mkdocs", *args)
