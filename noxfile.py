"""Automated testing linting and formatting apparatus."""

from __future__ import annotations

import nox
from nox.sessions import Session

package = "morfeus"
nox.options.sessions = "lint", "tests", "mypy", "deptry"  # default session
locations = "morfeus", "tests", "noxfile.py"  # Linting locations
pyversions = ["3.10"]


# Testing
@nox.session(python=pyversions)
def tests(session: Session) -> None:
    """Run tests."""
    args = session.posargs + ["--cov=morfeus", "--import-mode=importlib", "-s"]
    session.install("pytest", "pytest-cov")
    session.install(".")
    session.run("pytest", *args)


# Dependency checking
@nox.session(python=pyversions)
def deptry(session: Session) -> None:
    """Check for missing or unused dependencies."""
    session.install("deptry")
    session.install(".")
    session.run("deptry", "morfeus", "tests")


# Linting
@nox.session(python="3.10")
def lint(session: Session) -> None:
    """Lint code."""
    args = session.posargs or locations
    session.install(
        "flake8",
        "flake8-black",
        "flake8-bugbear",
        "flake8-import-order",
        "flake8-annotations",
        "flake8-docstrings",
        "darglint",
    )
    session.run("flake8", *args)


# Code formatting
@nox.session(python="3.10")
def black(session: Session) -> None:
    """Format code."""
    args = session.posargs or locations
    session.install("black")
    session.run("black", *args)


# Static typing
@nox.session(python="3.10")
def mypy(session: Session) -> None:
    """Run the static type checker."""
    args = session.posargs or locations
    session.install("mypy")
    session.install("numpy")
    session.install("types-setuptools")
    session.run("mypy", *args)
