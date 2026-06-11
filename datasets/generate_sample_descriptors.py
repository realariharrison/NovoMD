"""Generate the sample descriptor dataset.

Runs the open NovoMD calculator over a public list of well-known molecules and
writes outcome-level descriptors to ``sample_descriptors.csv``. Everything here
is public: the SMILES are textbook structures, and the descriptors are exactly
what ``pip install novomd`` produces. No proprietary inputs, models, or fields.

Reproduce:
    pip install novomd
    python datasets/generate_sample_descriptors.py
"""

from __future__ import annotations

import csv
import os
from typing import List, Tuple

from novomd import calculate_properties_batch

# (name, SMILES, category). Public, non-proprietary structures only.
COMPOUNDS: List[Tuple[str, str, str]] = [
    # Analgesics and anti-inflammatories
    ("aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O", "analgesic"),
    ("ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "analgesic"),
    ("acetaminophen", "CC(=O)NC1=CC=C(C=C1)O", "analgesic"),
    ("naproxen", "COC1=CC2=CC=C(C=C2C=C1)C(C)C(=O)O", "analgesic"),
    ("diclofenac", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "analgesic"),
    ("ketoprofen", "CC(C(=O)O)c1cccc(c1)C(=O)c1ccccc1", "analgesic"),
    # Stimulants and methylxanthines
    ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "stimulant"),
    ("theobromine", "CN1C=NC2=C1C(=O)NC(=O)N2C", "stimulant"),
    ("theophylline", "CN1C2=C(C(=O)N(C1=O)C)NC=N2", "stimulant"),
    ("nicotine", "CN1CCCC1C2=CN=CC=C2", "stimulant"),
    # Other common drugs
    ("metformin", "CN(C)C(=N)NC(=N)N", "antidiabetic"),
    ("warfarin", "CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O", "anticoagulant"),
    ("diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", "anxiolytic"),
    ("ibuprofen_lysine", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "analgesic"),
    ("salicylic_acid", "OC(=O)c1ccccc1O", "analgesic"),
    # Amino acids
    ("glycine", "C(C(=O)O)N", "amino_acid"),
    ("alanine", "CC(C(=O)O)N", "amino_acid"),
    ("serine", "C(C(C(=O)O)N)O", "amino_acid"),
    ("valine", "CC(C)C(C(=O)O)N", "amino_acid"),
    ("leucine", "CC(C)CC(C(=O)O)N", "amino_acid"),
    ("proline", "C1CC(NC1)C(=O)O", "amino_acid"),
    ("phenylalanine", "c1ccc(cc1)CC(C(=O)O)N", "amino_acid"),
    ("tyrosine", "c1cc(ccc1CC(C(=O)O)N)O", "amino_acid"),
    ("tryptophan", "c1ccc2c(c1)c(c[nH]2)CC(C(=O)O)N", "amino_acid"),
    ("cysteine", "C(C(C(=O)O)N)S", "amino_acid"),
    ("methionine", "CSCCC(C(=O)O)N", "amino_acid"),
    ("aspartic_acid", "C(C(C(=O)O)N)C(=O)O", "amino_acid"),
    ("glutamic_acid", "C(CC(=O)O)C(C(=O)O)N", "amino_acid"),
    ("lysine", "C(CCN)CC(C(=O)O)N", "amino_acid"),
    ("histidine", "c1c(nc[nH]1)CC(C(=O)O)N", "amino_acid"),
    ("threonine", "CC(C(C(=O)O)N)O", "amino_acid"),
    # Simple organics
    ("methane", "C", "organic"),
    ("ethanol", "CCO", "organic"),
    ("methanol", "CO", "organic"),
    ("benzene", "c1ccccc1", "organic"),
    ("toluene", "Cc1ccccc1", "organic"),
    ("phenol", "c1ccc(cc1)O", "organic"),
    ("acetic_acid", "CC(=O)O", "organic"),
    ("acetone", "CC(=O)C", "organic"),
    ("formaldehyde", "C=O", "organic"),
    ("urea", "C(=O)(N)N", "organic"),
    ("ethylene_glycol", "C(CO)O", "organic"),
    ("glycerol", "C(C(CO)O)O", "organic"),
    ("glucose", "C(C1C(C(C(C(O1)O)O)O)O)O", "sugar"),
    # Vitamins
    ("ascorbic_acid", "C(C(C1C(=C(C(=O)O1)O)O)O)O", "vitamin"),
    ("niacin", "O=C(O)c1cccnc1", "vitamin"),
    ("pyridoxine", "Cc1ncc(c(c1O)CO)CO", "vitamin"),
    # Neurotransmitters
    ("dopamine", "c1cc(c(cc1CCN)O)O", "neurotransmitter"),
    ("serotonin", "c1cc2c(cc1O)c(c[nH]2)CCN", "neurotransmitter"),
    ("histamine", "c1c(nc[nH]1)CCN", "neurotransmitter"),
    ("gaba", "C(CC(=O)O)CN", "neurotransmitter"),
    ("adrenaline", "CNCC(c1ccc(c(c1)O)O)O", "neurotransmitter"),
]

# Per-atom list fields are omitted from the flat table.
_LIST_FIELDS = {"coords_x", "coords_y", "coords_z", "atom_types", "bonds"}

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, "sample_descriptors.csv")


def main() -> None:
    smiles = [s for _, s, _ in COMPOUNDS]
    results = calculate_properties_batch(smiles)

    rows = []
    for (name, smi, category), result in zip(COMPOUNDS, results):
        if result["status"] != "ok":
            print(f"skipped {name}: {result['error']}")
            continue
        row = {"name": name, "category": category, "smiles": smi}
        for key, value in result["properties"].items():
            if key not in _LIST_FIELDS and key != "smiles":
                row[key] = value
        rows.append(row)

    fieldnames = list(rows[0].keys())
    with open(OUTPUT, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
