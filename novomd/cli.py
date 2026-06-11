"""Command-line entry point for NovoMD.

Currently exposes ``novomd props`` for local property calculation. Batch input
and CSV output (``novomd batch``) land in a follow-up change.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import List, Optional

from .__about__ import __version__
from .exceptions import NovoMDError


def _cmd_props(args: argparse.Namespace) -> int:
    from .core import calculate_properties

    try:
        result = calculate_properties(
            args.smiles,
            add_hydrogens=not args.no_hydrogens,
            optimize_3d=not args.no_optimize,
        )
    except NovoMDError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    indent = None if args.compact else 2
    print(json.dumps(result, indent=indent))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="novomd",
        description="Local-first molecular property calculator.",
    )
    parser.add_argument("--version", action="version", version=f"novomd {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    props = subparsers.add_parser("props", help="Compute descriptors for one SMILES string.")
    props.add_argument("smiles", help="SMILES string, e.g. 'CCO'")
    props.add_argument("--no-hydrogens", action="store_true", help="Do not add explicit hydrogens.")
    props.add_argument("--no-optimize", action="store_true", help="Skip 3D geometry optimization.")
    props.add_argument(
        "--compact", action="store_true", help="Emit single-line JSON instead of indented."
    )
    props.set_defaults(func=_cmd_props)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    exit_code: int = args.func(args)
    return exit_code


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
