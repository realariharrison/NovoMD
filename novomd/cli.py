"""Command-line entry point for NovoMD.

- ``novomd props "<smiles>"`` computes descriptors for one molecule.
- ``novomd batch <file.smi> --out results.csv`` processes many at once.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from typing import Any, Dict, List, Optional

from .__about__ import __version__
from .exceptions import NovoMDError

# Property fields that are per-atom lists; omitted from flat CSV output.
_LIST_FIELDS = {"coords_x", "coords_y", "coords_z", "atom_types", "bonds"}


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


def _read_smiles_file(path: str) -> List[str]:
    """Read a .smi file: one SMILES per line; the first whitespace-separated
    token is taken as the SMILES. Blank lines and ``#`` comments are skipped."""
    molecules: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            molecules.append(stripped.split()[0])
    return molecules


def _flatten_for_csv(results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for item in results:
        row: Dict[str, Any] = {
            "smiles": item["smiles"],
            "status": item["status"],
            "error": item.get("error", ""),
        }
        for key, value in item.get("properties", {}).items():
            if key not in _LIST_FIELDS:
                row[key] = value
        rows.append(row)
    return rows


def _cmd_batch(args: argparse.Namespace) -> int:
    from .batch import calculate_properties_batch

    try:
        molecules = _read_smiles_file(args.input)
    except OSError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if not molecules:
        print(f"error: no SMILES found in {args.input}", file=sys.stderr)
        return 1

    results = calculate_properties_batch(
        molecules,
        add_hydrogens=not args.no_hydrogens,
        optimize_3d=not args.no_optimize,
    )
    succeeded = sum(1 for r in results if r["status"] == "ok")
    print(
        f"processed {len(results)} molecules: {succeeded} ok, {len(results) - succeeded} failed",
        file=sys.stderr,
    )

    if args.out:
        if args.out.lower().endswith((".csv", ".tsv")):
            rows = _flatten_for_csv(results)
            fieldnames: List[str] = []
            for row in rows:
                for key in row:
                    if key not in fieldnames:
                        fieldnames.append(key)
            delimiter = "\t" if args.out.lower().endswith(".tsv") else ","
            with open(args.out, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
                writer.writeheader()
                writer.writerows(rows)
        else:
            with open(args.out, "w", encoding="utf-8") as handle:
                json.dump(results, handle, indent=2)
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        print(json.dumps(results, indent=2))

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

    batch = subparsers.add_parser(
        "batch", help="Compute descriptors for many SMILES from a .smi file."
    )
    batch.add_argument("input", help="Path to a .smi file (one SMILES per line).")
    batch.add_argument(
        "--out",
        help="Write results to this file (.csv/.tsv for a table, otherwise JSON). "
        "Without --out, JSON is printed to stdout.",
    )
    batch.add_argument("--no-hydrogens", action="store_true", help="Do not add explicit hydrogens.")
    batch.add_argument("--no-optimize", action="store_true", help="Skip 3D geometry optimization.")
    batch.set_defaults(func=_cmd_batch)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    exit_code: int = args.func(args)
    return exit_code


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
