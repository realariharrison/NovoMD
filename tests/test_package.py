"""Tests for the importable `novomd` library (local-first, no server)."""

import pytest


def test_version_exposed():
    import novomd

    assert isinstance(novomd.__version__, str)
    assert novomd.__version__.count(".") >= 2


def test_public_api_surface():
    import novomd

    for name in ("calculate_properties", "smiles_to_pdb", "pdb_to_omd", "get_atom_type"):
        assert hasattr(novomd, name)


def test_invalid_smiles_raises():
    from novomd import InvalidSMILESError, calculate_properties

    with pytest.raises(InvalidSMILESError):
        calculate_properties("INVALID_SMILES_STRING")


@pytest.mark.skipif(
    not pytest.importorskip("rdkit", reason="RDKit not installed"),
    reason="RDKit required for property calculation",
)
class TestCalculateProperties:
    """End-to-end local calculation: the P0 acceptance path."""

    def test_ethanol_returns_descriptor_dict(self):
        from novomd import calculate_properties

        result = calculate_properties("CCO")

        assert isinstance(result, dict)
        assert result["smiles"] == "CCO"
        # CCO = C2H6O, MW ~46.07
        assert 45.0 < result["molecular_weight"] < 47.0
        assert result["num_heavy_atoms"] == 3  # two carbons + one oxygen
        # descriptor families are all present
        for key in (
            "radius_of_gyration",
            "conformer_energy",
            "dipole_moment",
            "molecular_volume",
            "atom_types",
        ):
            assert key in result

    def test_no_hydrogens_flag(self):
        from novomd import calculate_properties

        result = calculate_properties("C", add_hydrogens=False)
        assert result["num_heavy_atoms"] == 1
