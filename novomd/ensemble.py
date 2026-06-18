"""Optional conformer-ensemble generation via openconf.

This module is the *only* place NovoMD touches openconf. openconf is an optional
dependency (``pip install 'novomd[ensemble]'``); when it is absent, or when the
running Python is older than openconf supports (it requires 3.12+), the public
entry points here degrade gracefully so the base, single-conformer calculator
keeps working unchanged.

openconf (MIT, https://github.com/rowansci/openconf) generates the geometries;
NovoMD computes its descriptors on top of them. We use openconf's own MMFF94s
energies and Boltzmann weights rather than re-deriving them.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np

try:
    from openconf import generate_conformers, preset_config  # type: ignore

    OPENCONF_AVAILABLE = True
except Exception:  # pragma: no cover - exercised only without openconf / on py<3.12
    OPENCONF_AVAILABLE = False


# openconf presets that make sense for property prediction. "ensemble" is the
# openconf-recommended default for this use case; "macrocycle" enables low-mode
# following plus a wide energy window for large rings.
_ALLOWED_PRESETS = {"rapid", "ensemble", "spectroscopic", "docking", "analogue", "macrocycle"}


class EnsembleUnavailableError(RuntimeError):
    """Raised when an ensemble is requested but cannot be produced.

    In normal use this is caught and the caller falls back to the single
    conformer path. It only propagates when ``strict_ensemble=True``.
    """


@dataclass
class RawEnsemble:
    """Geometry-level ensemble handed back to the core descriptor routine."""

    coords: List[np.ndarray]  # per-conformer (N, 3) float arrays, Angstrom
    atoms: List[str]  # element symbols, shared across conformers
    pdb_blocks: List[str]  # per-conformer PDB block (for partial-charge parsing)
    weights: np.ndarray  # Boltzmann weights, length == n_conformers
    energies: List[float]  # openconf MMFF94s energies (kcal/mol)
    n_conformers: int


def ensemble_available() -> bool:
    """Return True if openconf can be imported in this environment."""
    return OPENCONF_AVAILABLE


def generate_ensemble(
    smiles: str,
    *,
    max_conformers: int = 50,
    preset: str = "ensemble",
    random_seed: int = 42,
    temperature: float = 298.15,
) -> RawEnsemble:
    """Generate a conformer ensemble for a SMILES string with openconf.

    Args:
        smiles: molecule as a SMILES string.
        max_conformers: cap on returned conformers (openconf ``max_out``).
        preset: one of ``rapid``, ``ensemble``, ``spectroscopic``, ``docking``,
            ``analogue``, ``macrocycle``.
        random_seed: reproducibility seed.
        temperature: temperature (K) for Boltzmann weighting.

    Raises:
        EnsembleUnavailableError: openconf is not importable or produced nothing.
        InvalidSMILESError: the SMILES could not be parsed.
        ValueError: unknown preset.
    """
    if not OPENCONF_AVAILABLE:
        raise EnsembleUnavailableError(
            "openconf is not installed. Install the ensemble extra with: "
            "pip install 'novomd[ensemble]' (requires Python 3.12+)."
        )

    # Imported here so this module loads cleanly even without RDKit present.
    from rdkit import Chem

    from .exceptions import InvalidSMILESError

    if preset not in _ALLOWED_PRESETS:
        raise ValueError(
            f"Unknown openconf preset {preset!r}; choose one of {sorted(_ALLOWED_PRESETS)}."
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSMILESError(f"Invalid SMILES string: {smiles!r}")

    # openconf forbids passing both a config and a preset, so start from the
    # preset config and override the fields we care about.
    config = preset_config(preset)
    config.max_out = max_conformers
    config.random_seed = random_seed

    ens = generate_conformers(mol, config=config)

    n = ens.n_conformers
    if n == 0:
        raise EnsembleUnavailableError(f"openconf produced no conformers for {smiles!r}.")

    atoms = [a.GetSymbol() for a in ens.mol.GetAtoms()]
    conf_ids = ens.conf_ids
    coords = [np.asarray(ens.coords(i), dtype=float) for i in range(n)]
    pdb_blocks = [str(Chem.MolToPDBBlock(ens.mol, confId=conf_ids[i]) or "") for i in range(n)]
    weights = np.asarray(ens.boltzmann_weights(temperature=temperature), dtype=float)
    energies = list(ens.energies)

    return RawEnsemble(
        coords=coords,
        atoms=atoms,
        pdb_blocks=pdb_blocks,
        weights=weights,
        energies=energies,
        n_conformers=n,
    )
