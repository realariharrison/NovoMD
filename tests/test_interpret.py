"""Tests for the interpretation layer (drug-likeness + summary)."""

import pytest

rdkit = pytest.importorskip("rdkit", reason="RDKit required for interpretation")


class TestCalculateDruglikeness:
    def test_keys_present(self):
        from novomd import calculate_druglikeness

        d = calculate_druglikeness("CCO")
        for key in (
            "molecular_weight",
            "logp",
            "tpsa",
            "h_bond_donors",
            "h_bond_acceptors",
            "rotatable_bonds",
            "aromatic_rings",
            "fraction_csp3",
            "qed",
            "lipinski",
            "veber",
        ):
            assert key in d

    def test_ethanol_is_drug_like(self):
        from novomd import calculate_druglikeness

        d = calculate_druglikeness("CCO")
        assert 45 < d["molecular_weight"] < 47
        assert d["lipinski"]["violations"] == []
        assert d["lipinski"]["within_ro5"] is True
        assert d["veber"]["passes"] is True
        assert 0.0 <= d["qed"] <= 1.0

    def test_large_lipophilic_molecule_flags_violations(self):
        from novomd import calculate_druglikeness

        # A long-chain fatty acid: very lipophilic and floppy.
        d = calculate_druglikeness("CCCCCCCCCCCCCCCCCCCC(=O)O")
        assert "logP > 5" in d["lipinski"]["violations"]
        assert d["veber"]["passes"] is False

    def test_invalid_smiles_raises(self):
        from novomd import InvalidSMILESError, calculate_druglikeness

        with pytest.raises(InvalidSMILESError):
            calculate_druglikeness("NOT_A_SMILES")


class TestSummarize:
    def test_summary_is_descriptive_string(self):
        from novomd import calculate_druglikeness, summarize

        text = summarize(calculate_druglikeness("CC(=O)OC1=CC=CC=C1C(=O)O"))
        assert isinstance(text, str)
        assert "Lipinski" in text
        assert "QED" in text


class TestInterpret:
    def test_interpret_combines_data_and_summary(self):
        from novomd import interpret

        result = interpret("CCO")
        assert "summary" in result
        assert "molecular_weight" in result
        assert isinstance(result["summary"], str)


class TestExplainCLI:
    def test_explain_table(self, capsys):
        from novomd.cli import main

        rc = main(["explain", "CCO"])
        out = capsys.readouterr().out
        assert rc == 0
        assert "QED" in out
        assert "Lipinski" in out
        assert "molecular weight" in out

    def test_explain_json(self, capsys):
        import json

        from novomd.cli import main

        rc = main(["explain", "CCO", "--json"])
        out = capsys.readouterr().out
        assert rc == 0
        payload = json.loads(out)
        assert payload["smiles"] == "CCO"
        assert "summary" in payload
