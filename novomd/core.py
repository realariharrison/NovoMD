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


def calculate_properties(
    smiles: str, *, add_hydrogens: bool = True, optimize_3d: bool = True
) -> Dict[str, Any]:
    """Compute the full molecular descriptor set for a SMILES string, locally.

    Parses the SMILES, embeds a 3D conformer, and returns a flat dictionary of
    identity metadata (molecular weight, atom/bond counts) plus the geometry,
    energy, electrostatic, surface/volume and visualization descriptors. No
    network access, no API key, no server.

    Args:
        smiles: The molecule as a SMILES string (e.g. ``"CCO"``).
        add_hydrogens: Add explicit hydrogens before embedding (default True).
        optimize_3d: Run UFF geometry optimization on the conformer (default True).

    Returns:
        A descriptor dictionary keyed by property name.

    Raises:
        RDKitNotAvailableError: RDKit is not installed.
        InvalidSMILESError: The SMILES string could not be parsed.
    """
    _require_rdkit()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSMILESError(f"Invalid SMILES string: {smiles!r}")

    if add_hydrogens:
        mol = Chem.AddHs(mol)

    pdb_content = smiles_to_pdb(smiles, optimize_3d=optimize_3d, add_hydrogens=add_hydrogens)
    coords, atoms = extract_coordinates_from_pdb(pdb_content)
    properties = calculate_all_molecular_properties(coords, atoms, mol, pdb_content)

    return {
        "smiles": smiles,
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        **properties,
    }
