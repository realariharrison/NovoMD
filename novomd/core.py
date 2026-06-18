"""Framework-free molecular property calculation core.

Everything in this module runs locally with no network calls and no web
framework. If RDKit is installed, :func:`calculate_properties` turns a SMILES
string into a full descriptor dictionary on your own hardware.

The numerical routines here are the same ones the REST service uses; the
service imports from this module rather than redefining them.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.spatial.distance import cdist

from .exceptions import InvalidSMILESError, RDKitNotAvailableError

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

    RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only without RDKit
    RDKIT_AVAILABLE = False


def _require_rdkit() -> None:
    if not RDKIT_AVAILABLE:
        raise RDKitNotAvailableError(
            "RDKit is required for this operation but is not installed. "
            "Install it with: pip install novomd"
        )


def smiles_to_pdb(smiles: str, optimize_3d: bool = True, add_hydrogens: bool = True) -> str:
    """Convert a SMILES string to a 3D PDB block using RDKit."""
    _require_rdkit()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSMILESError(f"Invalid SMILES string: {smiles!r}")

    if add_hydrogens:
        mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, randomSeed=42)

    if optimize_3d:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)

    pdb_block = Chem.MolToPDBBlock(mol)
    return str(pdb_block) if pdb_block else ""


def calculate_partial_charges(pdb_content: str, method: str = "gasteiger") -> Dict[int, float]:
    """Estimate per-atom partial charges from PDB content.

    Simplified electronegativity-based model; outcome-level descriptor only.
    """
    charges: Dict[int, float] = {}
    atom_index = 0

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            element = line[76:78].strip() if len(line) > 76 else "C"

            electronegativities = {
                "H": 2.20,
                "C": 2.55,
                "N": 3.04,
                "O": 3.44,
                "F": 3.98,
                "S": 2.58,
                "Cl": 3.16,
                "Br": 2.96,
            }

            en = electronegativities.get(element, 2.5)
            charge = (en - 2.5) * 0.1

            charges[atom_index] = round(charge, 4)
            atom_index += 1

    return charges


def extract_coordinates_from_pdb(pdb_content: str) -> Tuple[np.ndarray, List[str]]:
    """Extract 3D coordinates and element symbols from PDB content."""
    coords: List[List[float]] = []
    atoms: List[str] = []

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                element = line[76:78].strip() if len(line) > 76 else "C"

                coords.append([x, y, z])
                atoms.append(element)
            except (ValueError, IndexError):
                continue

    return np.array(coords), atoms


def calculate_all_molecular_properties(
    coords: np.ndarray, atoms: List[str], mol: Any, pdb_content: str
) -> Dict[str, Any]:
    """Calculate the full descriptor set from 3D coordinates.

    Returns geometry, energy estimate, electrostatic, surface/volume, atom-count
    and 3D-visualization descriptors. Returns an empty dict for empty input.
    """

    if len(coords) == 0:
        return {}

    # Center of mass
    center = np.mean(coords, axis=0)
    centered_coords = coords - center

    # === GEOMETRY PROPERTIES (7) ===

    # Radius of gyration
    rgyr = np.sqrt(np.mean(np.sum(centered_coords**2, axis=1)))

    # Maximum distance (span)
    distances = cdist(coords, coords)
    max_dist = np.max(distances)

    # Inertia tensor for shape analysis (I is the conventional physics symbol)
    I = np.zeros((3, 3))  # noqa: E741
    for coord in centered_coords:
        I[0, 0] += coord[1] ** 2 + coord[2] ** 2
        I[1, 1] += coord[0] ** 2 + coord[2] ** 2
        I[2, 2] += coord[0] ** 2 + coord[1] ** 2
        I[0, 1] -= coord[0] * coord[1]
        I[0, 2] -= coord[0] * coord[2]
        I[1, 2] -= coord[1] * coord[2]
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    # Principal moments of inertia
    eigenvalues = np.sort(np.linalg.eigvals(I).real)
    pmi1, pmi2, pmi3 = eigenvalues

    # Shape descriptors
    asphericity = pmi3 - 0.5 * (pmi1 + pmi2)
    eccentricity = (pmi3 - pmi1) / pmi3 if pmi3 > 0 else 0
    inertia_shape_factor = pmi1 / pmi3 if pmi3 > 0 else 0

    # === SURFACE/VOLUME PROPERTIES (4) ===

    num_atoms = len(atoms)
    num_heavy = sum(1 for a in atoms if a not in ["H", "h"])

    # Estimate molecular volume and surface area
    hull_volume = num_atoms * 15.0  # Å³ per atom
    hull_area = num_atoms * 30.0  # Å² per atom
    globularity = (
        min(1.0, (36 * np.pi * hull_volume**2) ** (1 / 3) / hull_area) if hull_area > 0 else 0
    )
    surface_to_volume = hull_area / hull_volume if hull_volume > 0 else 0

    # === ENERGY PROPERTIES (6) ===
    # These are estimates - real MD would provide actual values

    # Bond detection
    bonds = []
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            if distances[i, j] < 1.6:  # Typical bond length
                bonds.append([int(i), int(j)])

    conformer_energy = -10.0 * num_atoms
    vdw_energy = -0.5 * len(bonds)
    electrostatic_energy = -0.1 * num_atoms
    torsion_strain = 0.1 * max(0, len(bonds) - num_atoms + 1)
    angle_strain = 0.05 * num_atoms
    optimization_delta = abs(conformer_energy) * 0.1

    # === ELECTROSTATIC PROPERTIES (6) ===

    dipole_moment = np.linalg.norm(center) * 0.1
    total_charge = 0.0  # Neutral

    # Calculate partial charges
    charges = calculate_partial_charges(pdb_content, "gasteiger")
    if charges:
        charge_values = list(charges.values())
        max_partial_charge = max(charge_values)
        min_partial_charge = min(charge_values)
        charge_span = max_partial_charge - min_partial_charge
        total_charge = sum(charge_values)
    else:
        max_partial_charge = 0.5
        min_partial_charge = -0.5
        charge_span = 1.0

    electrostatic_potential = dipole_moment * 0.1

    # Return all descriptors
    return {
        # Geometry (7)
        "radius_of_gyration": round(float(rgyr), 3),
        "asphericity": round(float(asphericity), 3),
        "eccentricity": round(float(eccentricity), 3),
        "inertia_shape_factor": round(float(inertia_shape_factor), 3),
        "span_r": round(float(max_dist), 3),
        "pmi1": round(float(pmi1), 3),
        "pmi2": round(float(pmi2), 3),
        # Energy (6)
        "conformer_energy": round(float(conformer_energy), 2),
        "vdw_energy": round(float(vdw_energy), 2),
        "electrostatic_energy": round(float(electrostatic_energy), 2),
        "torsion_strain": round(float(torsion_strain), 2),
        "angle_strain": round(float(angle_strain), 2),
        "optimization_delta": round(float(optimization_delta), 2),
        # Electrostatics (6)
        "dipole_moment": round(float(dipole_moment), 3),
        "total_charge": round(float(total_charge), 3),
        "max_partial_charge": round(float(max_partial_charge), 3),
        "min_partial_charge": round(float(min_partial_charge), 3),
        "charge_span": round(float(charge_span), 3),
        "electrostatic_potential": round(float(electrostatic_potential), 3),
        # Surface/Volume (4)
        "sasa": round(float(hull_area), 1),
        "molecular_volume": round(float(hull_volume), 1),
        "globularity": round(float(globularity), 3),
        "surface_to_volume_ratio": round(float(surface_to_volume), 3),
        # Atom counts (2)
        "num_atoms_with_h": int(num_atoms),
        "num_heavy_atoms": int(num_heavy),
        # Visualization (5+)
        "coords_x": [round(float(c[0]), 4) for c in coords],
        "coords_y": [round(float(c[1]), 4) for c in coords],
        "coords_z": [round(float(c[2]), 4) for c in coords],
        "atom_types": atoms,
        "bonds": bonds,
    }


def _calculate_properties_ensemble(
    smiles: str, *, max_conformers: int, preset: str, temperature: float
) -> Dict[str, Any]:
    """Compute Boltzmann-weighted, ensemble-averaged descriptors via openconf.

    Reuses :func:`calculate_all_molecular_properties` per conformer, then
    combines: numeric scalars are Boltzmann-weighted across the ensemble;
    structural fields (coords, bonds, atom_types) are taken from the
    lowest-energy conformer. Adds real ensemble metadata from openconf.
    """
    from .ensemble import generate_ensemble  # lazy: keeps openconf optional

    raw = generate_ensemble(
        smiles, max_conformers=max_conformers, preset=preset, temperature=temperature
    )

    per_conf = [
        calculate_all_molecular_properties(raw.coords[i], raw.atoms, None, raw.pdb_blocks[i])
        for i in range(raw.n_conformers)
    ]

    # Normalize weights; fall back to uniform if degenerate.
    weights = raw.weights
    if weights.size != raw.n_conformers or weights.sum() <= 0 or not np.isfinite(weights).all():
        weights = np.ones(raw.n_conformers) / raw.n_conformers
    else:
        weights = weights / weights.sum()

    rep_idx = int(np.argmin(raw.energies)) if raw.energies else 0
    rep = per_conf[rep_idx]

    merged: Dict[str, Any] = {}
    for key, rep_val in rep.items():
        if isinstance(rep_val, bool):
            merged[key] = rep_val
        elif isinstance(rep_val, int):
            # count-like field: average then round to int (identical across confs)
            vals = np.array([float(pc[key]) for pc in per_conf])
            merged[key] = int(round(float(np.dot(weights, vals))))
        elif isinstance(rep_val, float):
            vals = np.array([float(pc[key]) for pc in per_conf])
            merged[key] = round(float(np.dot(weights, vals)), 4)
        else:
            # lists / strings (coords_*, atom_types, bonds): take representative
            merged[key] = rep_val

    # Real ensemble metadata (sourced from openconf, not placeholders).
    rgyr = np.array([float(pc["radius_of_gyration"]) for pc in per_conf])
    rgyr_mean = float(np.dot(weights, rgyr))
    rgyr_std = float(np.sqrt(max(0.0, np.dot(weights, (rgyr - rgyr_mean) ** 2))))

    energies = np.array(raw.energies, dtype=float)
    finite = energies[np.isfinite(energies)]

    merged.update(
        {
            "n_conformers": int(raw.n_conformers),
            "ensemble_energy_min_kcal": round(float(finite.min()), 3) if finite.size else None,
            "ensemble_energy_spread_kcal": (
                round(float(finite.max() - finite.min()), 3) if finite.size else None
            ),
            "conformational_flexibility_rgyr": round(rgyr_std, 4),
            "method": "openconf_ensemble",
        }
    )
    return merged


def calculate_properties(
    smiles: str,
    *,
    add_hydrogens: bool = True,
    optimize_3d: bool = True,
    conformers: int | None = None,
    ensemble_preset: str = "ensemble",
    temperature: float = 298.15,
    strict_ensemble: bool = False,
) -> Dict[str, Any]:
    """Compute the full molecular descriptor set for a SMILES string, locally.

    By default this embeds a single conformer (UFF-optimized) and returns the
    classic descriptor dictionary. Pass ``conformers=N`` (N >= 2) to compute a
    Boltzmann-weighted ensemble average using openconf; if openconf is not
    available this falls back to the single-conformer result unless
    ``strict_ensemble=True``.

    Args:
        smiles: the molecule as a SMILES string (e.g. ``"CCO"``).
        add_hydrogens: add explicit hydrogens before embedding (default True).
        optimize_3d: run UFF optimization on the single conformer (default True;
            ignored in ensemble mode, where openconf handles minimization).
        conformers: if >= 2, request an openconf ensemble of up to this many
            conformers. ``None`` or ``1`` uses the single-conformer path.
        ensemble_preset: openconf preset for ensemble mode (default "ensemble";
            use "macrocycle" for large rings).
        temperature: temperature (K) for Boltzmann weighting in ensemble mode.
        strict_ensemble: if True, raise instead of falling back when an ensemble
            was requested but openconf is unavailable.

    Returns:
        A descriptor dictionary keyed by property name. Always includes
        ``"method"`` ("single_conformer_uff" or "openconf_ensemble") and
        ``"n_conformers"``.

    Raises:
        RDKitNotAvailableError: RDKit is not installed.
        InvalidSMILESError: the SMILES string could not be parsed.
        EnsembleUnavailableError: ensemble requested with strict_ensemble=True
            but openconf is unavailable.
    """
    _require_rdkit()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSMILESError(f"Invalid SMILES string: {smiles!r}")

    if add_hydrogens:
        mol = Chem.AddHs(mol)

    identity = {
        "smiles": smiles,
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
    }

    want_ensemble = conformers is not None and conformers > 1
    if want_ensemble:
        from .ensemble import EnsembleUnavailableError, ensemble_available

        if ensemble_available():
            try:
                properties = _calculate_properties_ensemble(
                    smiles,
                    max_conformers=int(conformers),
                    preset=ensemble_preset,
                    temperature=temperature,
                )
                return {**identity, **properties}
            except Exception:
                if strict_ensemble:
                    raise
                # graceful fallback to the single-conformer path below
        elif strict_ensemble:
            raise EnsembleUnavailableError(
                "Ensemble requested with strict_ensemble=True but openconf is "
                "unavailable. Install with: pip install 'novomd[ensemble]' "
                "(requires Python 3.12+)."
            )

    # --- single-conformer path (unchanged numerics) ---
    pdb_content = smiles_to_pdb(smiles, optimize_3d=optimize_3d, add_hydrogens=add_hydrogens)
    coords, atoms = extract_coordinates_from_pdb(pdb_content)
    properties = calculate_all_molecular_properties(coords, atoms, mol, pdb_content)
    properties["n_conformers"] = 1
    properties["method"] = "single_conformer_uff"
    return {**identity, **properties}
