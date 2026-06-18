"""Batch molecular property calculation with per-item error isolation.

Process a list of SMILES in one call. A single malformed or un-embeddable
molecule never fails the whole batch: each item returns its own status, so a
1,000-molecule run with a few bad entries still returns every good result.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

from .core import calculate_properties

# Default ceiling on how many molecules one call will process. The REST service
# enforces the same cap so a single request cannot tie up the worker.
MAX_BATCH_SIZE = 1000


def calculate_properties_batch(
    smiles_list: Sequence[str],
    *,
    add_hydrogens: bool = True,
    optimize_3d: bool = True,
    conformers: Optional[int] = None,
    ensemble_preset: str = "ensemble",
    temperature: float = 298.15,
    strict_ensemble: bool = False,
    max_batch_size: Optional[int] = MAX_BATCH_SIZE,
) -> List[Dict[str, Any]]:
    """Compute descriptors for many SMILES, isolating per-item failures.

    Args:
        smiles_list: SMILES strings to process.
        add_hydrogens: Add explicit hydrogens before embedding (default True).
        optimize_3d: Run UFF geometry optimization on each conformer (default True).
        max_batch_size: Reject inputs larger than this. Pass ``None`` to disable
            the check (the caller is then responsible for bounding the work).

    Returns:
        One result dict per input, in order. Each is either
        ``{"smiles": ..., "status": "ok", "properties": {...}}`` or
        ``{"smiles": ..., "status": "error", "error": "<message>"}``.

    Raises:
        ValueError: The input is larger than ``max_batch_size``.
    """
    if max_batch_size is not None and len(smiles_list) > max_batch_size:
        raise ValueError(f"Batch size {len(smiles_list)} exceeds the maximum of {max_batch_size}.")

    results: List[Dict[str, Any]] = []
    for smiles in smiles_list:
        try:
            properties = calculate_properties(
                smiles,
                add_hydrogens=add_hydrogens,
                optimize_3d=optimize_3d,
                conformers=conformers,
                ensemble_preset=ensemble_preset,
                temperature=temperature,
                strict_ensemble=strict_ensemble,
            )
            results.append({"smiles": smiles, "status": "ok", "properties": properties})
        except Exception as exc:  # noqa: BLE001 - one bad molecule must not kill the batch
            results.append({"smiles": smiles, "status": "error", "error": str(exc)})

    return results
