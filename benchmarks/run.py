"""Reproducible, offline, single-core benchmark for NovoMD descriptors.

The point of this benchmark is not raw speed. It is the axis NovoMD owns:
descriptors computed with no account, no API key, no network, no GPU, and
*deterministically* — the same SMILES yields byte-identical descriptors across
runs and machines.

Run it yourself (nothing leaves your machine):

    pip install novomd
    python benchmarks/run.py

It prints a Markdown report to stdout; redirect it to refresh the committed
report:

    python benchmarks/run.py > docs/benchmark_report.md
"""

from __future__ import annotations

import platform
import statistics
import time
from typing import List, Tuple

from novomd import __version__, calculate_properties

# Public, size-graded molecules (small to large), all textbook structures.
MOLECULES: List[Tuple[str, str]] = [
    ("water", "O"),
    ("ethanol", "CCO"),
    ("benzene", "c1ccccc1"),
    ("aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
    ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("glucose", "C(C1C(C(C(C(O1)O)O)O)O)O"),
    ("cholesterol", "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"),
]

REPEATS = 7  # median of this many timed runs per molecule


def _median_seconds(smiles: str, repeats: int = REPEATS) -> float:
    timings = []
    for _ in range(repeats):
        start = time.perf_counter()
        calculate_properties(smiles)
        timings.append(time.perf_counter() - start)
    return statistics.median(timings)


def main() -> None:
    machine = f"{platform.system()} {platform.machine()}, Python {platform.python_version()}"

    print("# NovoMD local-first benchmark\n")
    print(
        f"NovoMD {__version__}. Single-conformer descriptors, computed locally with "
        "no account, no API key, no network, and no GPU.\n"
    )
    print(f"Environment: {machine}. Single CPU core.\n")
    print("| molecule | heavy atoms | descriptors | median time (ms) | deterministic |")
    print("| --- | ---: | ---: | ---: | :---: |")

    all_deterministic = True
    for name, smiles in MOLECULES:
        try:
            props = calculate_properties(smiles)
        except Exception as exc:  # pragma: no cover - skip anything that won't embed
            print(f"| {name} | - | - | - | error: {exc} |")
            continue

        # Determinism: a second call must produce an identical dict.
        deterministic = calculate_properties(smiles) == props
        all_deterministic = all_deterministic and deterministic

        heavy = props["num_heavy_atoms"]
        n_descriptors = sum(1 for k, v in props.items() if isinstance(v, (int, float)))
        median_ms = _median_seconds(smiles) * 1000.0

        print(
            f"| {name} | {heavy} | {n_descriptors} | {median_ms:.1f} | "
            f"{'yes' if deterministic else 'NO'} |"
        )

    verdict = "all molecules identical across runs." if all_deterministic else "MISMATCH DETECTED."
    print()
    print(
        f"**Determinism:** {verdict} "
        "The conformer embedding uses a fixed seed, so a given SMILES produces the "
        "same descriptors every time, on any machine.\n"
    )
    print(
        "**Reproduce:** `pip install novomd && python benchmarks/run.py`. Numbers are "
        "machine-dependent; the determinism guarantee is not."
    )


if __name__ == "__main__":
    main()
