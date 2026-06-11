# Contributing to NovoMD

Thank you for your interest in contributing to NovoMD! This document provides guidelines and instructions for contributing.

## Code of Conduct

Be respectful, inclusive, and professional in all interactions.

## Scientific Integrity Guidelines

### No Mock Data Policy

NovoMD is a **scientific computing tool** used for real computational chemistry research. To maintain trust and scientific validity:

**STRICTLY PROHIBITED:**
- ❌ Mock data generation functions
- ❌ Random number generators for scientific properties
- ❌ Simulated/fake calculation results
- ❌ Placeholder data in API responses

**REQUIRED:**
- ✅ All molecular properties must be calculated from real 3D structures
- ✅ Use established scientific libraries (RDKit, numpy, scipy)
- ✅ Document calculation methods and assumptions
- ✅ Cite algorithms and force field parameters used
- ✅ If a calculation isn't implemented, return null/error - never fake it

### Adding New Endpoints

When proposing new endpoints:

1. **Provide Scientific Justification**
   - What property/calculation does it provide?
   - What is the scientific use case?
   - What algorithms/methods will be used?

2. **Show Real Implementation**
   - No "TODO: implement actual calculation" with mock data
   - If the full implementation requires external tools (GROMACS, AMBER), document integration requirements
   - Stub endpoints must return errors, not simulated data

3. **Include Validation**
   - How will results be validated?
   - What are known test cases?
   - Compare with established tools when possible

### Example: Good vs Bad

❌ **BAD - Mock Data:**
```python
@app.post("/calculate-logp")
def calculate_logp(smiles: str):
    # TODO: implement real calculation
    return {"logp": random.uniform(0, 5)}
```

✅ **GOOD - Real Calculation:**
```python
@app.post("/calculate-logp")
def calculate_logp(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(400, "Invalid SMILES")
    logp = Descriptors.MolLogP(mol)  # RDKit's real calculation
    return {"logp": round(logp, 2), "method": "Wildman-Crippen"}
```

✅ **ACCEPTABLE - Not Yet Implemented:**
```python
@app.post("/calculate-binding-affinity")
def calculate_binding_affinity(protein: str, ligand: str):
    raise HTTPException(
        501,
        "Binding affinity requires MD simulation engine integration. "
        "Please integrate with GROMACS or AMBER for production use."
    )
```

### Property Calculation Standards

All molecular properties must:
- Be derived from actual molecular structures (3D coordinates, topology)
- Use peer-reviewed algorithms or established software packages
- Include units in response
- Document limitations and assumptions

## How to Contribute

### Reporting Bugs

1. Check if the bug has already been reported in [Issues](https://github.com/realariharrison/NovoMD/issues)
2. If not, create a new issue with:
   - Clear, descriptive title
   - Steps to reproduce
   - Expected vs actual behavior
   - Environment details (OS, Python version, Docker version)
   - Relevant logs or error messages

### Suggesting Features

1. Check if the feature has been requested in [Issues](https://github.com/realariharrison/NovoMD/issues)
2. Create a new issue with:
   - Clear description of the feature
   - Use case and motivation
   - Proposed implementation (if applicable)

### Pull Requests

1. **Fork the repository**
   ```bash
   git clone https://github.com/YOUR_USERNAME/NovoMD.git
   cd NovoMD
   ```

2. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Make your changes**
   - Follow the code style guidelines below
   - Add tests for new functionality
   - Update documentation as needed

4. **Test your changes**
   ```bash
   # Run the service locally
   uvicorn main:app --reload

   # Test endpoints
   python test_novomd.py
   ```

5. **Commit your changes**
   ```bash
   git add .
   git commit -m "Add: brief description of your changes"
   ```

   Commit message format:
   - `Add:` for new features
   - `Fix:` for bug fixes
   - `Update:` for improvements to existing features
   - `Docs:` for documentation changes
   - `Refactor:` for code refactoring

6. **Push and create a pull request**
   ```bash
   git push origin feature/your-feature-name
   ```

   Then create a PR on GitHub with:
   - Clear description of changes
   - Link to related issues
   - Screenshots (if applicable)

## Code Style Guidelines

### Python

- Follow [PEP 8](https://pep8.org/) style guide
- Use type hints for function signatures
- Write docstrings for all public functions
- Keep functions focused and concise
- Use meaningful variable names

Example:
```python
def calculate_binding_energy(
    protein_pdb: str,
    ligand_smiles: str,
    method: str = "MM-GBSA"
) -> float:
    """
    Calculate protein-ligand binding energy.

    Args:
        protein_pdb: PDB ID or content of the protein
        ligand_smiles: SMILES string of the ligand
        method: Calculation method (MM-GBSA, MM-PBSA, etc.)

    Returns:
        Binding energy in kcal/mol

    Raises:
        ValueError: If invalid inputs provided
    """
    # Implementation here
    pass
```

### API Endpoints

- Use clear, RESTful endpoint names
- Include comprehensive docstrings
- Define Pydantic models for request/response schemas
- Add appropriate error handling
- Use proper HTTP status codes

### Documentation

- Keep README.md up to date
- Document all API endpoints
- Include usage examples
- Update CHANGELOG.md for significant changes

## Development Setup

1. **Install development dependencies**
   ```bash
   pip install -r requirements.txt
   pip install -r requirements-dev.txt  # If available
   ```

2. **Set up pre-commit hooks** (recommended)
   ```bash
   pip install pre-commit
   pre-commit install
   ```

3. **Run in development mode**
   ```bash
   uvicorn main:app --reload --log-level debug
   ```

## Testing

- Write tests for all new features
- Ensure existing tests pass
- Aim for high code coverage
- Test edge cases and error conditions

```bash
# Run tests
pytest tests/

# Run with coverage
pytest --cov=. tests/
```

## Documentation

- Update API documentation in README.md
- Add inline code comments for complex logic
- Update docstrings when modifying functions
- Keep examples current and working

## Review Process

1. Automated checks must pass (linting, tests)
2. At least one maintainer review required
3. Address all review comments
4. Squash commits before merge (if requested)

## Release Process

Maintainers will handle releases following semantic versioning:
- **MAJOR**: Breaking changes
- **MINOR**: New features (backward compatible)
- **PATCH**: Bug fixes

## Questions?

- Open an issue for general questions
- Use [Discussions](https://github.com/realariharrison/NovoMD/discussions) for broader topics
- Tag maintainers if urgent: @maintainer-username

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to NovoMD!
