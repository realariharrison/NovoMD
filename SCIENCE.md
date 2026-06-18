# How NovoMD computes its descriptors

This document describes the methods behind NovoMD, and is honest about what is
real physics and what is a deterministic estimate. NovoMD is a descriptor
calculator, not a simulation engine; the value it offers is being local,
deterministic, and clear about its own boundary.

## Single-conformer method (default)

For each SMILES, NovoMD:

1. Parses the molecule with RDKit and adds explicit hydrogens.
2. Embeds one 3D conformer with RDKit's distance-geometry embedder
   (`EmbedMolecule`, ETKDG defaults) using a **fixed random seed (42)**.
3. Optimizes that conformer with the UFF force field.
4. Computes descriptors from the resulting coordinates.

What is computed from the real 3D geometry:

- **Geometry / shape**: radius of gyration, asphericity, eccentricity, inertia
  shape factor, span, principal moments of inertia.
- **Surface / volume**: SASA, molecular volume, globularity, surface-to-volume
  ratio (atom-count estimates scaled from the geometry).
- **Electrostatics**: dipole estimate and a simple electronegativity-based
  partial-charge spread.
- **Atom / bond counts** and the full coordinate, atom-type, and bond lists.

**Honest caveat.** The fields named `conformer_energy`, `vdw_energy`,
`electrostatic_energy`, `torsion_strain`, and `angle_strain` in the
single-conformer output are **deterministic estimates derived from atom and
bond counts, not force-field energies**. They are useful as stable relative
features; they are not physical energies. For real conformer energies, use the
ensemble path below (or a dedicated engine).

## Drug-likeness (topological)

`calculate_druglikeness` / `novomd explain` compute standard, public
cheminformatics directly from the molecular graph (no 3D needed): molecular
weight, logP, TPSA, H-bond donor/acceptor counts, rotatable bonds, aromatic
rings, fraction sp3, QED, and the Lipinski (rule of five) and Veber verdicts.
These are rules of thumb, reported descriptively.

## Optional ensemble method (`novomd[ensemble]`)

When you pass `conformers=N`, NovoMD delegates **geometry generation** to
[openconf](https://github.com/rowansci/openconf) (MIT, by the Rowan team), then
computes its own descriptors on the resulting population:

- **Geometries and energies** come from openconf (MMFF94s).
- **Conformer-dependent scalars** (shape, span, dipole, the per-conformer
  descriptor set) are **Boltzmann-weighted** across the ensemble at the
  requested temperature, using openconf's weights.
- **Structural fields** (coordinates, atom types, bonds) are taken from the
  **lowest-energy** conformer, so they remain a consistent single structure.
- **Real ensemble energetics** are surfaced as new, clearly-named fields:
  `ensemble_energy_min_kcal`, `ensemble_energy_spread_kcal`, and
  `conformational_flexibility_rgyr` (the Boltzmann-weighted spread of the radius
  of gyration). These come from openconf, not from the estimates above.

We use openconf for geometry and Boltzmann weighting, not for ranking by MMFF
score: openconf's own documentation notes that MMFF rankings are weak, so we
lean on the weights and geometries, which is the supported use. NovoMD does not
re-derive conformers or dress up its estimate fields as ensemble physics.

## Determinism and reproducibility

The single-conformer embedding uses a fixed seed, so a given SMILES produces
**byte-identical descriptors every time, on any machine**, with no network call
and no account. The ensemble path is likewise seeded and deterministic for a
fixed openconf version. See [`docs/benchmark_report.md`](docs/benchmark_report.md)
for a reproducible, single-core, offline benchmark, and `benchmarks/run.py` to
run it yourself.

## What NovoMD is not

Deliberately. NovoMD describes a molecule; it does not predict its behavior.

- It does **not** run molecular dynamics trajectories.
- It does **not** dock, score binding affinity, or run FEP.
- It does **not** predict ADMET, pKa, solubility, or toxicity.

That boundary is the design. For predictive and simulation work, the same team
builds NovoMCP.
