"""Tests for batch property calculation (library, REST endpoint, CLI)."""

import pytest

rdkit = pytest.importorskip("rdkit", reason="RDKit required for property calculation")


class TestCalculatePropertiesBatch:
    """Library-level batch with per-item error isolation."""

    def test_all_valid(self):
        from novomd import calculate_properties_batch

        results = calculate_properties_batch(["CCO", "C", "O"])

        assert len(results) == 3
        assert all(r["status"] == "ok" for r in results)
        assert results[0]["smiles"] == "CCO"
        assert "properties" in results[0]

    def test_one_bad_smiles_isolated(self):
        """A malformed SMILES must not fail the whole batch."""
        from novomd import calculate_properties_batch

        results = calculate_properties_batch(["CCO", "NOT_A_SMILES", "C"])

        assert len(results) == 3
        assert [r["status"] for r in results] == ["ok", "error", "ok"]
        assert "error" in results[1]
        assert "properties" not in results[1]

    def test_hundred_molecules_one_malformed(self):
        """Brief acceptance: 100 molecules, one bad -> 99 ok + 1 error."""
        from novomd import calculate_properties_batch

        molecules = ["CCO"] * 99 + ["THIS_IS_NOT_VALID"]
        results = calculate_properties_batch(molecules)

        assert len(results) == 100
        assert sum(1 for r in results if r["status"] == "ok") == 99
        assert sum(1 for r in results if r["status"] == "error") == 1

    def test_empty_list(self):
        from novomd import calculate_properties_batch

        assert calculate_properties_batch([]) == []

    def test_exceeds_cap_raises(self):
        from novomd import calculate_properties_batch

        with pytest.raises(ValueError, match="exceeds the maximum"):
            calculate_properties_batch(["C", "C", "C"], max_batch_size=2)

    def test_cap_disabled(self):
        from novomd import calculate_properties_batch

        results = calculate_properties_batch(["C", "C", "C"], max_batch_size=None)
        assert len(results) == 3


class TestBatchEndpoint:
    """Tests for POST /batch."""

    def test_batch_requires_auth(self, client):
        response = client.post("/batch", json={"molecules": ["CCO"]})
        assert response.status_code == 401

    def test_empty_molecules_returns_400(self, client, auth_headers):
        response = client.post("/batch", json={"molecules": []}, headers=auth_headers)
        assert response.status_code == 400

    def test_oversized_batch_returns_400(self, client, auth_headers):
        molecules = ["CCO"] * 1001
        response = client.post("/batch", json={"molecules": molecules}, headers=auth_headers)
        assert response.status_code == 400
        assert "exceeds" in response.json()["detail"]

    def test_batch_with_one_malformed_returns_200(self, client, auth_headers):
        """99 valid + 1 malformed -> HTTP 200, 99 succeeded, 1 failed."""
        molecules = ["CCO"] * 99 + ["NOT_VALID_SMILES"]
        response = client.post("/batch", json={"molecules": molecules}, headers=auth_headers)

        assert response.status_code == 200
        data = response.json()
        assert data["count"] == 100
        assert data["succeeded"] == 99
        assert data["failed"] == 1
        assert len(data["results"]) == 100
        # the malformed entry carries a structured error, not properties
        bad = next(r for r in data["results"] if r["status"] == "error")
        assert "error" in bad and "properties" not in bad


class TestBatchCLI:
    """Tests for `novomd batch`."""

    def test_batch_writes_csv(self, tmp_path):
        from novomd.cli import main

        smi = tmp_path / "molecules.smi"
        smi.write_text("CCO\nC ethane_like\n# a comment\nNOT_VALID\n")
        out = tmp_path / "results.csv"

        exit_code = main(["batch", str(smi), "--out", str(out)])

        assert exit_code == 0
        text = out.read_text()
        # header + 3 data rows (comment line skipped)
        lines = [ln for ln in text.splitlines() if ln.strip()]
        assert lines[0].startswith("smiles,status,error")
        assert len(lines) == 1 + 3
        assert "molecular_weight" in lines[0]
        assert "NOT_VALID,error" in text
