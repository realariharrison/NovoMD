"""PDB to OpenMD (.omd) format conversion.

Framework-free helpers for turning a PDB block into an OpenMD input file and
mapping elements onto force-field atom types. Used by both the library and the
REST service.
"""

from __future__ import annotations


def get_atom_type(element: str, force_field: str) -> str:
    """Map an element symbol to a force-field atom type."""

    force_field_mappings = {
        "AMBER": {
            "H": "HC",
            "C": "CT",
            "N": "N",
            "O": "O",
            "S": "S",
            "P": "P",
            "F": "F",
            "Cl": "Cl",
            "Br": "Br",
        },
        "CHARMM": {
            "H": "HGA1",
            "C": "CG321",
            "N": "NG321",
            "O": "OG311",
            "S": "SG311",
            "P": "PG1",
            "F": "FGA1",
            "Cl": "CLGA1",
            "Br": "BRGA1",
        },
        "OPLS": {
            "H": "opls_140",
            "C": "opls_135",
            "N": "opls_238",
            "O": "opls_236",
            "S": "opls_200",
            "P": "opls_393",
            "F": "opls_164",
            "Cl": "opls_151",
            "Br": "opls_156",
        },
    }

    mapping = force_field_mappings.get(force_field, force_field_mappings["AMBER"])
    return mapping.get(element, element)


def pdb_to_omd(pdb_content: str, force_field: str, box_size: float, charge_method: str) -> str:
    """Convert a PDB block to OpenMD (.omd) format."""

    # Parse PDB content to extract atoms
    atoms = []
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_info = {
                "index": int(line[6:11].strip()),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "x": float(line[30:38].strip()),
                "y": float(line[38:46].strip()),
                "z": float(line[46:54].strip()),
                "element": line[76:78].strip() if len(line) > 76 else "C",
            }
            atoms.append(atom_info)

    if not atoms:
        raise ValueError("No atoms found in PDB content")

    # Generate OpenMD format content
    omd_content = """<OpenMD version=2>
  <MetaData>
    <molecule id="0">
      <name>Converted_Molecule</name>"""

    # Add atom definitions
    for atom in atoms:
        # Assign atom type based on element and force field
        atom_type = get_atom_type(str(atom["element"]), force_field)
        omd_content += f"""
      <atom id="{atom['index']}">
        <type>{atom_type}</type>
        <position x="{atom['x']}" y="{atom['y']}" z="{atom['z']}"/>
      </atom>"""

    omd_content += f"""
    </molecule>

    <forceField>{force_field}</forceField>
    <ensemble>NVT</ensemble>
    <target_temp>300</target_temp>
    <target_pressure>1</target_pressure>
  </MetaData>

  <Snapshot>
    <FrameData>
      <Time>0</Time>
      <Hmat>
        <Hxx>{box_size}</Hxx>
        <Hxy>0</Hxy>
        <Hxz>0</Hxz>
        <Hyx>0</Hyx>
        <Hyy>{box_size}</Hyy>
        <Hyz>0</Hyz>
        <Hzx>0</Hzx>
        <Hzy>0</Hzy>
        <Hzz>{box_size}</Hzz>
      </Hmat>
    </FrameData>

    <StuntDoubles>"""

    # Add positions for each atom
    for atom in atoms:
        atom_index: int = atom["index"]  # type: ignore[assignment]
        omd_content += f"""
      <StuntDouble index="{atom_index - 1}">
        <position x="{atom['x']}" y="{atom['y']}" z="{atom['z']}"/>
        <velocity x="0" y="0" z="0"/>
      </StuntDouble>"""

    omd_content += """
    </StuntDoubles>
  </Snapshot>
</OpenMD>"""

    return omd_content
