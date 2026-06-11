# Changelog

All notable changes to NovoMD will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.3.0] - 2026-06-11

### Added
- Interpretation layer: `calculate_druglikeness()` computes the standard
  medicinal-chemistry descriptors (logP, TPSA, H-bond donors/acceptors,
  rotatable bonds, aromatic rings, fraction sp3, QED) and the Lipinski and
  Veber rule-of-thumb checks.
- `summarize()` renders a plain-language profile; `interpret()` returns both.
- `novomd explain "<smiles>"` prints the drug-likeness table and summary
  (`--json` for machine output).

## [1.2.1] - 2026-06-11

### Added
- Designed CLI presentation: a bare `novomd` now shows a branded panel (with a
  terracotta accent on a terminal, plain text when piped) instead of an
  argparse error.

## [1.2.0] - 2026-06-11

### Added
- Importable `novomd` library: `calculate_properties()` computes the full
  descriptor set locally from a SMILES string, with no server and no API key.
- Batch processing with per-item error isolation: `calculate_properties_batch()`
  and a `POST /batch` endpoint (capped at 1,000 molecules per call).
- `novomd` command-line interface (`novomd props`, `novomd batch`).
- PyPI packaging with a `[server]` extra; published via Trusted Publishing.
- Sample descriptor dataset under `datasets/` as an open give-back.

### Changed
- The FastAPI service now imports the shared framework-free core rather than
  redefining the calculation logic.
- README rewritten around the local-first library, with a brand release card.

## [1.1.0] - 2025-01-08

### Added
- Comprehensive test suite with pytest (70%+ coverage target)
- CI/CD pipeline with GitHub Actions
  - Automated testing on Python 3.10, 3.11, 3.12
  - Code linting (flake8, black, isort)
  - Type checking (mypy)
  - Security scanning (bandit, safety)
  - Docker build verification
- `pyproject.toml` for modern Python packaging
- `requirements-dev.txt` for development dependencies
- Rate limiting middleware (100 requests/minute default)
- Pre-commit hooks configuration

### Changed
- Improved CORS configuration with environment variable control
- API key comparison now uses timing-attack resistant comparison
- Updated `/status` endpoint to reflect actual implemented capabilities

### Fixed
- Documentation inconsistencies between SECURITY.md and CONTRIBUTING.md
- Removed references to unimplemented features in API responses

### Security
- Added `secrets.compare_digest()` for API key validation
- Configurable CORS origins (no longer allows all origins by default)
- Added rate limiting to prevent abuse

## [1.0.0] - 2024-01-15

### Added
- Initial open-source release
- SMILES to OpenMD format conversion with 3D structure generation
- 32+ molecular property calculations from real 3D coordinates
  - Geometry properties (7): radius of gyration, asphericity, eccentricity, inertia shape factor, span, PMI1, PMI2
  - Energy properties (6): conformer energy, VDW energy, electrostatic energy, torsion strain, angle strain, optimization delta
  - Electrostatics properties (6): dipole moment, total charge, max/min partial charges, charge span, electrostatic potential
  - Surface/Volume properties (4): SASA, molecular volume, globularity, surface-to-volume ratio
  - Atom counts (2): total atoms, heavy atoms
  - 3D Visualization data (5+): coordinates, atom types, bonds
- Support for multiple force fields: AMBER, CHARMM, OPLS, GROMOS
- PDB to OpenMD conversion (`/atom2md` endpoint)
- API key authentication
- Docker and Docker Compose support
- Comprehensive documentation (README, CONTRIBUTING, SECURITY, DEPLOYMENT_GUIDE)

### Dependencies
- FastAPI 0.104.1
- Pydantic 2.5.0
- NumPy 1.24.3
- SciPy 1.11.4
- BioPython 1.81
- Optional: RDKit, OpenBabel

[Unreleased]: https://github.com/realariharrison/NovoMD/compare/v1.3.0...HEAD
[1.3.0]: https://github.com/realariharrison/NovoMD/compare/v1.2.1...v1.3.0
[1.2.1]: https://github.com/realariharrison/NovoMD/compare/v1.2.0...v1.2.1
[1.2.0]: https://github.com/realariharrison/NovoMD/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/realariharrison/NovoMD/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/realariharrison/NovoMD/releases/tag/v1.0.0
