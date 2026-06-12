"""
NovoMD Gradio interface and MCP server for Hugging Face Spaces.

A thin layer over the open-source ``novomd`` package: the same calculations that
``pip install novomd`` performs, exposed as a web UI and as MCP tools so AI
assistants can query them directly. The package is the single source of truth,
so the UI, the CLI, PyPI, and these MCP tools all return identical results.
"""

import json

import gradio as gr

from novomd import __version__, calculate_properties, generate_report, interpret
from novomd.exceptions import NovoMDError

try:
    from rdkit import Chem
    from rdkit.Chem import Draw

    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Force field options (affects the OpenMD output only, not property values)
FORCE_FIELDS = [
    ("AMBER14 - Recommended for proteins", "amber14"),
    ("AMBER99SB - Well-tested protein force field", "amber99sb"),
    ("CHARMM36 - Excellent for lipids", "charmm36"),
    ("OPLS-AA/M - Optimized for small molecules", "opls"),
    ("GROMOS 54A7 - United atom force field", "gromos54a7"),
]

# Example molecules
EXAMPLES = [
    ["CCO", "amber14"],  # Ethanol
    ["CC(=O)OC1=CC=CC=C1C(=O)O", "amber14"],  # Aspirin
    ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "amber14"],  # Caffeine
    ["CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "opls"],  # Ibuprofen
    ["c1ccccc1", "opls"],  # Benzene
]


def _properties_markdown(p: dict, force_field: str) -> str:
    return f"""
## Molecular Properties

### Basic Information
| Property | Value |
|----------|-------|
| SMILES | `{p['smiles']}` |
| Molecular Weight | {p['molecular_weight']:.2f} Da |
| Total Atoms (with H) | {p['num_atoms_with_h']} |
| Heavy Atoms | {p['num_heavy_atoms']} |
| Force Field | {force_field} |

### Geometry
| Property | Value |
|----------|-------|
| Radius of Gyration | {p['radius_of_gyration']:.3f} Å |
| Asphericity | {p['asphericity']:.3f} |
| Eccentricity | {p['eccentricity']:.3f} |
| Span (max distance) | {p['span_r']:.3f} Å |

### Surface & Volume
| Property | Value |
|----------|-------|
| SASA | {p['sasa']:.2f} Å² |
| Molecular Volume | {p['molecular_volume']:.2f} Å³ |
| Globularity | {p['globularity']:.3f} |
| Surface/Volume Ratio | {p['surface_to_volume_ratio']:.3f} |

### Energy (estimates)
| Property | Value |
|----------|-------|
| Conformer Energy | {p['conformer_energy']:.2f} |
| VDW Energy | {p['vdw_energy']:.2f} |
| Electrostatic Energy | {p['electrostatic_energy']:.2f} |
| Torsion Strain | {p['torsion_strain']:.2f} |
| Angle Strain | {p['angle_strain']:.2f} |

### Electrostatics
| Property | Value |
|----------|-------|
| Dipole Moment | {p['dipole_moment']:.3f} |
| Total Charge | {p['total_charge']:.4f} |
| Max Partial Charge | {p['max_partial_charge']:.4f} |
| Min Partial Charge | {p['min_partial_charge']:.4f} |
| Charge Span | {p['charge_span']:.4f} |

*Full 3D coordinates available in the JSON output below.*
"""


def process_molecule(smiles: str, force_field: str):
    """UI handler: returns a 2D image, a properties table, JSON, and any error."""
    if not smiles or not smiles.strip():
        return None, "Please enter a valid SMILES string.", "", ""

    smiles = smiles.strip()

    if not RDKIT_AVAILABLE:
        return None, "RDKit is not available. Please check the installation.", "", ""

    try:
        props = calculate_properties(smiles)
    except NovoMDError as exc:
        return None, f"Error processing molecule: {exc}", "", str(exc)

    mol_2d = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol_2d, size=(400, 300)) if mol_2d else None

    properties_md = _properties_markdown(props, force_field)
    json_output = json.dumps({"success": True, "force_field": force_field, **props}, indent=2)
    return img, properties_md, json_output, ""


# ---------------------------------------------------------------------------
# MCP tools — exposed to AI assistants. Each delegates to the novomd package,
# so the tools return exactly what the library and CLI return.
# ---------------------------------------------------------------------------


def molecular_properties(smiles: str) -> dict:
    """Compute 32+ 3D molecular descriptors for a SMILES string.

    Returns geometry (radius of gyration, asphericity, eccentricity, span,
    principal moments), an energy estimate, electrostatics (dipole, partial
    charges), surface and volume metrics, atom and bond counts, and the 3D
    coordinates. Computed locally with the open-source novomd package.

    Args:
        smiles: The molecule as a SMILES string, e.g. "CCO".
    """
    return calculate_properties(smiles)


def drug_likeness(smiles: str) -> dict:
    """Assess the drug-likeness of a SMILES string.

    Returns molecular weight, logP, TPSA, H-bond donor/acceptor counts,
    rotatable bonds, aromatic rings, fraction sp3, QED, the Lipinski (rule of
    five) and Veber verdicts, and a plain-language summary.

    This describes a molecule using public cheminformatics. It does not predict
    ADMET, pKa, solubility, or binding. For predictive work, see NovoMCP
    (novomcp.com).

    Args:
        smiles: The molecule as a SMILES string, e.g. "CCO".
    """
    return interpret(smiles)


def molecular_report(smiles: str, output_format: str = "markdown") -> str:
    """Generate a one-page molecular report for a SMILES string.

    Combines identity, drug-likeness, and a summary.

    Args:
        smiles: The molecule as a SMILES string, e.g. "CCO".
        output_format: "markdown" (default), "html" (with a 2D depiction), or "json".
    """
    return generate_report(smiles, fmt=output_format)


# Create Gradio interface
with gr.Blocks(title="NovoMD - Molecular Property Calculator") as demo:

    gr.Markdown(f"""
    # NovoMD - Molecular Property Calculator

    Calculate 32+ molecular properties and drug-likeness from SMILES strings,
    powered by the open-source [`novomd`](https://pypi.org/project/novomd/)
    package (v{__version__}).

    This describes molecules with public cheminformatics. It does not predict
    ADMET, pKa, solubility, or binding. For that, see NovoMCP.

    [GitHub](https://github.com/realariharrison/NovoMD) | [PyPI](https://pypi.org/project/novomd/)
    """)

    with gr.Row():
        with gr.Column(scale=1):
            smiles_input = gr.Textbox(
                label="SMILES String",
                placeholder="Enter a SMILES string (e.g., CCO for ethanol)",
                lines=1,
            )

            force_field_dropdown = gr.Dropdown(
                choices=FORCE_FIELDS,
                value="amber14",
                label="Force Field",
            )

            submit_btn = gr.Button("Calculate Properties", variant="primary")

            gr.Markdown("### Examples")
            gr.Examples(
                examples=EXAMPLES,
                inputs=[smiles_input, force_field_dropdown],
                label="Click to try:",
            )

            molecule_image = gr.Image(
                label="Molecule Structure",
                type="pil",
                elem_classes=["molecule-image"],
            )

        with gr.Column(scale=2):
            properties_output = gr.Markdown(
                label="Molecular Properties",
                value="Enter a SMILES string and click 'Calculate Properties' to see results.",
            )

            with gr.Accordion("JSON Output (for developers)", open=False):
                json_output = gr.Code(
                    label="JSON Response",
                    language="json",
                )

            error_output = gr.Textbox(
                label="Errors",
                visible=False,
            )

    # UI events. api_visibility="private" keeps these image-returning handlers out of the
    # MCP/API surface; the clean agent tools are registered via gr.api below.
    submit_btn.click(
        fn=process_molecule,
        inputs=[smiles_input, force_field_dropdown],
        outputs=[molecule_image, properties_output, json_output, error_output],
        api_visibility="private",
    )

    smiles_input.submit(
        fn=process_molecule,
        inputs=[smiles_input, force_field_dropdown],
        outputs=[molecule_image, properties_output, json_output, error_output],
        api_visibility="private",
    )

    # MCP tools (also callable as REST endpoints).
    gr.api(molecular_properties, api_name="molecular_properties")
    gr.api(drug_likeness, api_name="drug_likeness")
    gr.api(molecular_report, api_name="molecular_report")

    gr.Markdown("""
    ---

    **About NovoMD**

    NovoMD is an open-source, local-first molecular property calculator.
    Install it with `pip install novomd`, or run the REST service:

    ```bash
    docker run -d -p 8010:8010 -e NOVOMD_API_KEY="your-key" ghcr.io/realariharrison/novomd:latest
    ```

    MIT License | Built with RDKit and Gradio
    """)


if __name__ == "__main__":
    demo.launch(
        server_name="0.0.0.0",  # nosec B104
        server_port=7860,
        mcp_server=True,
    )
