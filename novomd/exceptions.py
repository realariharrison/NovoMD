"""Exception types raised by the NovoMD core.

These are plain Python exceptions with no web-framework coupling, so the core
can be imported and used as a library without FastAPI installed. The REST
service translates them into HTTP responses at the endpoint boundary.
"""

from __future__ import annotations


class NovoMDError(Exception):
    """Base class for all NovoMD errors."""


class InvalidSMILESError(NovoMDError, ValueError):
    """Raised when a SMILES string cannot be parsed into a molecule."""


class RDKitNotAvailableError(NovoMDError, RuntimeError):
    """Raised when an operation needs RDKit but it is not installed."""
