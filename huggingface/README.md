---
title: NovoMD
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: docker
app_port: 7860
pinned: false
license: mit
short_description: Molecular properties, drug-likeness, and reports from SMILES
tags:
  - chemistry
  - cheminformatics
  - rdkit
  - drug-discovery
  - mcp
  - computational-chemistry
---

# NovoMD - Molecular Property Calculator

Calculate molecular properties and drug-likeness from SMILES strings, powered by
the open-source [`novomd`](https://pypi.org/project/novomd/) package. The package
is the single source of truth, so this Space, the CLI, and the library all return
identical results.

## Features

- **Geometry**: radius of gyration, asphericity, eccentricity, span
- **Energy estimates**: conformer energy, VDW, electrostatic, strain
- **Electrostatics**: dipole moment, partial charges
- **Surface/volume**: SASA, molecular volume, globularity
- **Drug-likeness**: logP, TPSA, QED, Lipinski, Veber

This describes molecules with public cheminformatics. It does not predict ADMET,
pKa, solubility, or binding. For predictive work, see NovoMCP.

## MCP tools (for AI assistants)

This Space is also a Model Context Protocol server. Tools:

- `molecular_properties(smiles)` — 32+ 3D descriptors
- `drug_likeness(smiles)` — logP, TPSA, QED, Lipinski, Veber, plus a summary
- `molecular_report(smiles, output_format)` — a one-page report (markdown/html/json)

## Use the package directly

```bash
pip install novomd
```

## Links

- [PyPI](https://pypi.org/project/novomd/)
- [GitHub Repository](https://github.com/realariharrison/NovoMD)

## License

MIT License - Free for academic and commercial use.
