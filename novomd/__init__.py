"""NovoMD: a local-first molecular property calculator.

Compute molecular descriptors on your own hardware, no server and no API key::

    from novomd import calculate_properties
    props = calculate_properties("CCO")
    print(props["molecular_weight"])

The same core powers the optional REST service (``pip install novomd[server]``).
"""

from __future__ import annotations

from .__about__ import __version__
from .batch import MAX_BATCH_SIZE, calculate_properties_batch
from .conversion import get_atom_type, pdb_to_omd
from .core import (
    RDKIT_AVAILABLE,
    calculate_all_molecular_properties,
    calculate_partial_charges,
    calculate_properties,
    extract_coordinates_from_pdb,
    smiles_to_pdb,
)
from .exceptions import InvalidSMILESError, NovoMDError, RDKitNotAvailableError

__all__ = [
    "__version__",
    "RDKIT_AVAILABLE",
    "calculate_properties",
    "calculate_properties_batch",
    "MAX_BATCH_SIZE",
    "calculate_all_molecular_properties",
    "calculate_partial_charges",
    "extract_coordinates_from_pdb",
    "smiles_to_pdb",
    "get_atom_type",
    "pdb_to_omd",
    "NovoMDError",
    "InvalidSMILESError",
    "RDKitNotAvailableError",
]
