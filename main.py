import logging
import uuid
from datetime import datetime
from typing import Any, Dict, Optional

from fastapi import Depends, FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

# Rate limiting
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.errors import RateLimitExceeded
from slowapi.util import get_remote_address

# Molecular property core (framework-free, shared with the `novomd` library)
from novomd.conversion import get_atom_type, pdb_to_omd
from novomd.core import (
    RDKIT_AVAILABLE,
    calculate_all_molecular_properties,
    calculate_partial_charges,
    extract_coordinates_from_pdb,
    smiles_to_pdb,
)

# These core helpers are re-exported so existing callers and tests can keep
# importing them from `main` after the library extraction.
__all__ = [
    "app",
    "smiles_to_pdb",
    "pdb_to_omd",
    "get_atom_type",
    "calculate_partial_charges",
    "extract_coordinates_from_pdb",
    "calculate_all_molecular_properties",
]

# RDKit is optional at the service layer; the core handles its absence.
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:  # pragma: no cover
    Chem = None
    Descriptors = None

try:
    from openbabel import pybel  # noqa: F401  (imported to probe availability)

    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False
    logging.warning("OpenBabel not available - conversion features limited")

# Import config
from config import settings

# Configure logging
logging.basicConfig(
    level=getattr(logging, settings.log_level.upper()),
    format="[NovoMD] %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Initialize rate limiter
limiter = Limiter(key_func=get_remote_address)

# Initialize FastAPI
app = FastAPI(
    title="NovoMD - Molecular Dynamics Service",
    description="Open-source molecular dynamics simulation and docking service",
    version="1.1.0",
)

# Add rate limiter to app
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)  # type: ignore[arg-type]

# CORS middleware - use configurable origins
cors_origins = settings.get_cors_origins()
logger.info(f"CORS allowed origins: {cors_origins}")

app.add_middleware(
    CORSMiddleware,
    allow_origins=cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Authentication
from auth import verify_api_key


# Request/Response Models
class SMILESToOMDRequest(BaseModel):
    smiles: str = Field(..., description="SMILES string to convert")
    force_field: str = Field(default="AMBER", description="Force field: AMBER, CHARMM, OPLS")
    optimize_3d: bool = Field(default=True, description="Optimize 3D structure")
    add_hydrogens: bool = Field(default=True, description="Add explicit hydrogens")
    charge_method: str = Field(default="gasteiger", description="Charge calculation method")
    box_size: float = Field(default=30.0, description="Simulation box size in Angstroms")
    solvate: bool = Field(default=False, description="Add solvent molecules")


class OMDFileResponse(BaseModel):
    success: bool
    omd_content: Optional[str] = None
    pdb_content: Optional[str] = None
    metadata: Dict[str, Any]
    error: Optional[str] = None


# Helper Functions
def generate_job_id() -> str:
    """Generate unique job ID"""
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    unique_id = str(uuid.uuid4())[:6].upper()
    return f"MD_{timestamp}_{unique_id}"


# API Endpoints


@app.get("/health")
async def health(request: Request):
    """Health check endpoint"""
    return {
        "status": "healthy",
        "service": "NovoMD",
        "version": "1.0.0",
        "rdkit_available": RDKIT_AVAILABLE,
        "openbabel_available": OPENBABEL_AVAILABLE,
    }


@app.get("/status")
@limiter.limit(settings.rate_limit)
async def status(request: Request, api_key: str = Depends(verify_api_key)):
    """Service status and capabilities"""
    return {
        "service": "NovoMD",
        "version": "1.0.0",
        "capabilities": {
            "smiles_to_omd_conversion": True,
            "pdb_to_omd_conversion": True,
            "molecular_property_calculation": True,
            "geometry_analysis": True,
            "partial_charge_calculation": True,
            "rdkit_available": RDKIT_AVAILABLE,
            "openbabel_available": OPENBABEL_AVAILABLE,
        },
        "supported_force_fields": [
            "AMBER",
            "CHARMM",
            "OPLS",
        ],
        "property_categories": [
            "geometry",
            "energy",
            "electrostatics",
            "surface_volume",
            "atom_counts",
            "visualization",
        ],
    }


@app.get("/force-fields")
@limiter.limit(settings.rate_limit)
async def get_force_fields(request: Request, api_key: str = Depends(verify_api_key)):
    """Get available force fields and their descriptions"""
    return {
        "force_fields": [
            {
                "name": "amber14",
                "description": "AMBER ff14SB - recommended for proteins",
                "best_for": ["proteins", "nucleic acids"],
                "water_model": "TIP3P",
            },
            {
                "name": "amber99sb",
                "description": "AMBER ff99SB-ILDN - well-tested protein force field",
                "best_for": ["proteins"],
                "water_model": "TIP3P",
            },
            {
                "name": "charmm36",
                "description": "CHARMM36 - good for lipids and membranes",
                "best_for": ["lipids", "membranes", "proteins"],
                "water_model": "TIP3P",
            },
            {
                "name": "opls",
                "description": "OPLS-AA/M - optimized for small molecules",
                "best_for": ["small molecules", "organic compounds"],
                "water_model": "TIP4P",
            },
            {
                "name": "gromos54a7",
                "description": "GROMOS 54A7 - united atom force field",
                "best_for": ["proteins", "peptides"],
                "water_model": "SPC",
            },
        ]
    }


@app.post("/smiles-to-omd", response_model=OMDFileResponse)
@limiter.limit(settings.rate_limit)
async def convert_smiles_to_omd(
    request: Request, data: SMILESToOMDRequest, api_key: str = Depends(verify_api_key)
):
    """
    Convert SMILES string to OpenMD format (.omd file)

    This endpoint performs the complete conversion pipeline:
    1. SMILES → 3D structure (using RDKit)
    2. 3D structure → PDB format
    3. PDB → OpenMD format (simulating atom2md functionality)

    The generated .omd file includes:
    - Atomic coordinates
    - Force field parameters
    - Simulation box definition
    - Partial charges (if requested)
    """
    try:
        # Step 1: Convert SMILES to PDB
        logger.info(f"Converting SMILES to PDB: {data.smiles}")
        pdb_content = smiles_to_pdb(
            data.smiles, optimize_3d=data.optimize_3d, add_hydrogens=data.add_hydrogens
        )

        # Step 2: Convert PDB to OpenMD format
        logger.info(f"Converting PDB to OpenMD format with {data.force_field} force field")
        omd_content = pdb_to_omd(pdb_content, data.force_field, data.box_size, data.charge_method)

        # Step 3: Calculate all 32 molecular properties
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(data.smiles)
            if data.add_hydrogens:
                mol = Chem.AddHs(mol)

            # Extract coordinates from PDB
            coords, atoms = extract_coordinates_from_pdb(pdb_content)

            # Calculate comprehensive molecular properties
            properties = calculate_all_molecular_properties(coords, atoms, mol, pdb_content)

            # Build complete metadata with all 32+ properties
            metadata = {
                "smiles": data.smiles,
                "num_atoms": mol.GetNumAtoms(),
                "num_bonds": mol.GetNumBonds(),
                "molecular_weight": round(Descriptors.MolWt(mol), 2),
                "force_field": data.force_field,
                "box_size": data.box_size,
                "optimized": data.optimize_3d,
                "hydrogens_added": data.add_hydrogens,
                "charge_method": data.charge_method,
                "conversion_timestamp": datetime.now().isoformat(),
                **properties,  # Add all 32 calculated properties
            }
        else:
            metadata = {
                "smiles": data.smiles,
                "force_field": data.force_field,
                "box_size": data.box_size,
                "conversion_timestamp": datetime.now().isoformat(),
            }

        return OMDFileResponse(
            success=True, omd_content=omd_content, pdb_content=pdb_content, metadata=metadata
        )

    except HTTPException as e:
        raise e
    except Exception as e:
        logger.error(f"SMILES to OMD conversion failed: {str(e)}")
        return OMDFileResponse(success=False, error=str(e), metadata={"smiles": data.smiles})


class Atom2MDRequest(BaseModel):
    pdb_content: str = Field(..., description="PDB file content")
    force_field: str = Field(default="AMBER", description="Force field")
    box_size: float = Field(default=30.0, description="Box size in Angstroms")


@app.post("/atom2md")
@limiter.limit(settings.rate_limit)
async def atom2md_conversion(
    request: Request, data: Atom2MDRequest, api_key: str = Depends(verify_api_key)
):
    """
    Direct PDB to OpenMD conversion (atom2md equivalent)

    This endpoint simulates the atom2md tool functionality:
    - Takes PDB format input
    - Generates OpenMD format output
    - Assigns force field parameters
    - Sets up simulation box
    """
    try:
        omd_content = pdb_to_omd(data.pdb_content, data.force_field, data.box_size, "gasteiger")

        return {
            "success": True,
            "omd_content": omd_content,
            "metadata": {
                "force_field": data.force_field,
                "box_size": data.box_size,
                "conversion_timestamp": datetime.now().isoformat(),
            },
        }
    except Exception as e:
        logger.error(f"atom2md conversion failed: {str(e)}")
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/force-field-types/{force_field}")
@limiter.limit(settings.rate_limit)
async def get_force_field_atom_types(
    request: Request, force_field: str, api_key: str = Depends(verify_api_key)
):
    """Get available atom types for a specific force field"""
    force_field_types = {
        "AMBER": {
            "description": "Amber force field atom types",
            "common_types": {
                "HC": "H bonded to aliphatic carbon",
                "CT": "sp3 aliphatic carbon",
                "N": "sp2 nitrogen",
                "O": "sp2 oxygen",
                "OH": "hydroxyl oxygen",
                "S": "sulfur",
                "P": "phosphorus",
            },
        },
        "CHARMM": {
            "description": "CHARMM force field atom types",
            "common_types": {
                "HGA1": "aliphatic hydrogen",
                "CG321": "aliphatic carbon",
                "NG321": "neutral nitrogen",
                "OG311": "hydroxyl oxygen",
                "SG311": "sulfur",
            },
        },
        "OPLS": {
            "description": "OPLS-AA force field atom types",
            "common_types": {
                "opls_140": "alkane H",
                "opls_135": "alkane CH3",
                "opls_238": "amine nitrogen",
                "opls_236": "carbonyl oxygen",
            },
        },
    }

    if force_field not in force_field_types:
        raise HTTPException(status_code=404, detail=f"Force field {force_field} not found")

    return force_field_types[force_field]


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=settings.port)  # nosec B104
