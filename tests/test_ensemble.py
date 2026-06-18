"""Tests for the optional openconf ensemble path."""

from __future__ import annotations

import pytest

from novomd import calculate_properties
from novomd.ensemble import OPENCONF_AVAILABLE, EnsembleUnavailableError

pytest.importorskip("rdkit", reason="RDKit required for property calculation")

FLEXIBLE = "CCCCc1ccccc1OCCO"  # several rotatable bonds


def test_single_path_unchanged():
    """Default call still returns the classic shape plus the two new tags."""
    props = calculate_properties("CCO")
    assert props["method"] == "single_conformer_uff"
    assert props["n_conformers"] == 1
    assert props["molecular_weight"] == 46.07
    assert "radius_of_gyration" in props


def test_fallback_when_unavailable():
    """conformers>=2 must never crash; it falls back when openconf is absent."""
    props = calculate_properties(FLEXIBLE, conformers=20)
    if OPENCONF_AVAILABLE:
        assert props["method"] == "openconf_ensemble"
    else:
        assert props["method"] == "single_conformer_uff"
    # descriptor surface is intact either way
    assert "radius_of_gyration" in props


def test_strict_raises_when_unavailable():
    if OPENCONF_AVAILABLE:
        pytest.skip("openconf installed; strict mode would succeed")
    with pytest.raises(EnsembleUnavailableError):
        calculate_properties(FLEXIBLE, conformers=20, strict_ensemble=True)


@pytest.mark.skipif(not OPENCONF_AVAILABLE, reason="requires openconf")
def test_ensemble_fields_present():
    props = calculate_properties(FLEXIBLE, conformers=30)
    assert props["method"] == "openconf_ensemble"
    assert props["n_conformers"] >= 1
    assert "conformational_flexibility_rgyr" in props
    assert "ensemble_energy_spread_kcal" in props
    # ensemble output is a superset of the single-conformer keys
    single = calculate_properties(FLEXIBLE)
    assert set(single).issubset(set(props))


@pytest.mark.skipif(not OPENCONF_AVAILABLE, reason="requires openconf")
def test_batch_mixed_validity():
    from novomd import calculate_properties_batch

    results = calculate_properties_batch(["CCO", "NOT_VALID", FLEXIBLE], conformers=10)
    assert results[0]["status"] == "ok"
    assert results[1]["status"] == "error"
    assert results[2]["status"] == "ok"
