"""Command-line entry point for NovoMD.

- ``novomd props "<smiles>"`` computes descriptors for one molecule.
- ``novomd batch <file.smi> --out results.csv`` processes many at once.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from typing import Any, Dict, List, Optional, TextIO

from .__about__ import __version__
from .exceptions import NovoMDError

# Property fields that are per-atom lists; omitted from flat CSV output.
_LIST_FIELDS = {"coords_x", "coords_y", "coords_z", "atom_types", "bonds"}

# Brand terracotta (#B8704B) as a 24-bit ANSI color, used only on a TTY.
_TERRACOTTA = "\x1b[38;2;184;112;75m"
_BOLD = "\x1b[1m"
_DIM = "\x1b[2m"
_RESET = "\x1b[0m"


def _supports_color(stream: TextIO) -> bool:
    """Color only when writing to a real terminal and not opted out."""
    return (
        hasattr(stream, "isatty")
        and stream.isatty()
        and os.environ.get("NO_COLOR") is None
        and os.environ.get("TERM") != "dumb"
    )


def _render_panel(color: bool) -> str:
    """The designed presentation shown for a bare ``novomd`` invocation."""
    terra = _TERRACOTTA if color else ""
    bold = _BOLD if color else ""
    dim = _DIM if color else ""
    reset = _RESET if color else ""
    return f"""
 {terra}MOLECULAR PROPERTY CALCULATOR{reset}

 {terra}{bold}novomd{reset}  ·  v{__version__}

   {bold}props{reset}     descriptors for one molecule
             {dim}novomd props "CCO"{reset}

   {bold}explain{reset}   drug-likeness and a plain-language summary
             {dim}novomd explain "CCO"{reset}

   {bold}report{reset}    a one-page report (markdown, html, or json)
             {dim}novomd report "CCO" --out report.md{reset}

   {bold}batch{reset}     many molecules, one pass
             {dim}novomd batch mols.smi --out results.csv{reset}

 {dim}local-first · open source · novomcp.com{reset}
"""


def _cmd_props(args: argparse.Namespace) -> int:
    from .core import calculate_properties

    try:
        result = calculate_properties(
            args.smiles,
            add_hydrogens=not args.no_hydrogens,
            optimize_3d=not args.no_optimize,
            conformers=args.conformers,
            ensemble_preset=args.preset,
        )
    except NovoMDError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    # Note when an ensemble was requested but the run fell back to single-conformer.
    if args.conformers and args.conformers > 1 and result.get("method") != "openconf_ensemble":
        print(
            "note: openconf unavailable; computed a single conformer. "
            "Install with: pip install 'novomd[ensemble]' (Python 3.12+).",
            file=sys.stderr,
        )

    indent = None if args.compact else 2
    print(json.dumps(result, indent=indent))
    return 0


def _cmd_explain(args: argparse.Namespace) -> int:
    from .interpret import calculate_druglikeness, summarize

    try:
        data = calculate_druglikeness(args.smiles)
    except NovoMDError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps({**data, "summary": summarize(data)}, indent=2))
        return 0

    rows = [
        ("molecular weight", data["molecular_weight"]),
        ("logP", data["logp"]),
        ("TPSA", data["tpsa"]),
        ("H-bond donors", data["h_bond_donors"]),
        ("H-bond acceptors", data["h_bond_acceptors"]),
        ("rotatable bonds", data["rotatable_bonds"]),
        ("aromatic rings", data["aromatic_rings"]),
        ("fraction sp3", data["fraction_csp3"]),
        ("QED", data["qed"]),
    ]
    print(data["smiles"])
    print()
    for label, value in rows:
        print(f"  {label:<18} {value}")
    print()
    lipinski = data["lipinski"]["violations"]
    print(f"  Lipinski   {'no violations' if not lipinski else ', '.join(lipinski)}")
    veber = data["veber"]["violations"]
    print(f"  Veber      {'passes' if not veber else ', '.join(veber)}")
    print()
    print(summarize(data))
    return 0


def _cmd_report(args: argparse.Namespace) -> int:
    from .report import generate_report

    fmt = args.format
    if fmt is None:
        if args.out and args.out.lower().endswith(".html"):
            fmt = "html"
        elif args.out and args.out.lower().endswith(".json"):
            fmt = "json"
        else:
            fmt = "markdown"

    try:
        content = generate_report(args.smiles, fmt=fmt)
    except NovoMDError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            handle.write(content)
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        print(content)
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
        description="novomd · local-first molecular property calculation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"novomd {__version__}")

    # Not required: a bare `novomd` shows the panel instead of an argparse error.
    subparsers = parser.add_subparsers(dest="command", required=False, metavar="<command>")

    props = subparsers.add_parser("props", help="Compute descriptors for one SMILES string.")
    props.add_argument("smiles", help="SMILES string, e.g. 'CCO'")
    props.add_argument("--no-hydrogens", action="store_true", help="Do not add explicit hydrogens.")
    props.add_argument("--no-optimize", action="store_true", help="Skip 3D geometry optimization.")
    props.add_argument(
        "--conformers",
        type=int,
        default=None,
        metavar="N",
        help="Average over an openconf ensemble of up to N conformers "
        "(requires the [ensemble] extra; falls back to single conformer).",
    )
    props.add_argument(
        "--preset",
        default="ensemble",
        choices=["rapid", "ensemble", "spectroscopic", "docking", "analogue", "macrocycle"],
        help="openconf preset for ensemble mode (default: ensemble).",
    )
    props.add_argument(
        "--compact", action="store_true", help="Emit single-line JSON instead of indented."
    )
    props.set_defaults(func=_cmd_props)

    explain = subparsers.add_parser(
        "explain", help="Drug-likeness and a plain-language summary for one SMILES."
    )
    explain.add_argument("smiles", help="SMILES string, e.g. 'CCO'")
    explain.add_argument("--json", action="store_true", help="Emit JSON instead of a table.")
    explain.set_defaults(func=_cmd_explain)

    report = subparsers.add_parser(
        "report", help="Generate a one-page molecular report for one SMILES."
    )
    report.add_argument("smiles", help="SMILES string, e.g. 'CCO'")
    report.add_argument("--out", help="Write to this file (format inferred from .md/.html/.json).")
    report.add_argument(
        "--format",
        choices=["markdown", "html", "json"],
        help="Override the output format (default: inferred from --out, else markdown).",
    )
    report.set_defaults(func=_cmd_report)

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
    if getattr(args, "command", None) is None:
        # Bare `novomd`: show the designed panel rather than erroring.
        print(_render_panel(_supports_color(sys.stdout)))
        return 0
    exit_code: int = args.func(args)
    return exit_code


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
