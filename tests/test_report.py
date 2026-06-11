"""Tests for the one-page molecular report."""

import json

import pytest

rdkit = pytest.importorskip("rdkit", reason="RDKit required for reports")


class TestGenerateReport:
    def test_markdown(self):
        from novomd import generate_report

        md = generate_report("CC(=O)OC1=CC=CC=C1C(=O)O", "markdown")
        assert md.startswith("# Molecular report: C9H8O4")
        assert "| Molecular weight | 180.16 |" in md
        assert "Lipinski:" in md
        assert "QED" in md

    def test_json(self):
        from novomd import generate_report

        payload = json.loads(generate_report("CCO", "json"))
        assert payload["formula"] == "C2H6O"
        assert "summary" in payload
        assert "qed" in payload

    def test_html_has_depiction_and_brand(self):
        from novomd import generate_report

        html = generate_report("CCO", "html")
        assert "<svg" in html  # 2D structure depiction
        assert "C2H6O" in html
        assert "#B8704B" in html  # terracotta accent

    def test_invalid_format_raises(self):
        from novomd import generate_report

        with pytest.raises(ValueError, match="Unknown report format"):
            generate_report("CCO", "pdf")

    def test_invalid_smiles_raises(self):
        from novomd import InvalidSMILESError, generate_report

        with pytest.raises(InvalidSMILESError):
            generate_report("NOT_A_SMILES", "markdown")


class TestReportCLI:
    def test_writes_html_inferred_from_extension(self, tmp_path):
        from novomd.cli import main

        out = tmp_path / "aspirin.html"
        rc = main(["report", "CC(=O)OC1=CC=CC=C1C(=O)O", "--out", str(out)])
        assert rc == 0
        text = out.read_text()
        assert "<svg" in text and "C9H8O4" in text

    def test_format_override(self, tmp_path, capsys):
        from novomd.cli import main

        rc = main(["report", "CCO", "--format", "json"])
        assert rc == 0
        payload = json.loads(capsys.readouterr().out)
        assert payload["formula"] == "C2H6O"
