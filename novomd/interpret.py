"""Drug-likeness descriptors and plain-language interpretation.

Turns a SMILES string into the medicinal-chemistry descriptors people actually
screen on (logP, TPSA, H-bond counts, rotatable bonds, QED) and the textbook
rule-of-thumb checks (Lipinski's rule of five, Veber), then summarizes them in
a sentence.

These are standard, public cheminformatics computed with RDKit. This module
*describes* a molecule. It does not predict ADMET, pKa, solubility, or binding.
For predictive work, that boundary is where NovoMCP begins.
"""

from __future__ import annotations

from typing import Any, Dict

from .exceptions import InvalidSMILESError, RDKitNotAvailableError

try:
    from rdkit import Chem
    from rdkit.Chem import QED, Descriptors

    RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only without RDKit
    RDKIT_AVAILABLE = False


def calculate_druglikeness(smiles: str) -> Dict[str, Any]:
    """Compute medicinal-chemistry descriptors and rule-of-thumb verdicts.

    Returns a dict with molecular weight, logP, TPSA, H-bond donor/acceptor
    counts, rotatable bonds, aromatic rings, fraction sp3, QED, and the
    Lipinski and Veber assessments. Topological (2D) and fast: no 3D embedding.

    Raises:
        RDKitNotAvailableError: RDKit is not installed.
        InvalidSMILESError: The SMILES string could not be parsed.
    """
    if not RDKIT_AVAILABLE:
        raise RDKitNotAvailableError(
            "RDKit is required for this operation but is not installed. "
            "Install it with: pip install novomd"
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSMILESError(f"Invalid SMILES string: {smiles!r}")

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = int(Descriptors.NumHDonors(mol))
    hba = int(Descriptors.NumHAcceptors(mol))
    rotatable = int(Descriptors.NumRotatableBonds(mol))
    aromatic_rings = int(Descriptors.NumAromaticRings(mol))
    fraction_csp3 = Descriptors.FractionCSP3(mol)
    qed = QED.qed(mol)

    # Lipinski's rule of five. A molecule is conventionally "within" the rule
    # when it violates at most one of the four thresholds.
    lipinski_violations = []
    if mw > 500:
        lipinski_violations.append("MW > 500")
    if logp > 5:
        lipinski_violations.append("logP > 5")
    if hbd > 5:
        lipinski_violations.append("HBD > 5")
    if hba > 10:
        lipinski_violations.append("HBA > 10")

    # Veber criteria for oral bioavailability.
    veber_violations = []
    if rotatable > 10:
        veber_violations.append("rotatable bonds > 10")
    if tpsa > 140:
        veber_violations.append("TPSA > 140")

    return {
        "smiles": smiles,
        "molecular_weight": round(mw, 2),
        "logp": round(logp, 2),
        "tpsa": round(tpsa, 2),
        "h_bond_donors": hbd,
        "h_bond_acceptors": hba,
        "rotatable_bonds": rotatable,
        "aromatic_rings": aromatic_rings,
        "fraction_csp3": round(fraction_csp3, 3),
        "qed": round(qed, 3),
        "lipinski": {
            "violations": lipinski_violations,
            "within_ro5": len(lipinski_violations) <= 1,
        },
        "veber": {
            "violations": veber_violations,
            "passes": len(veber_violations) == 0,
        },
    }


def _size_word(mw: float) -> str:
    if mw < 150:
        return "very small"
    if mw < 350:
        return "small"
    if mw < 500:
        return "medium-sized"
    return "large"


def _lipophilicity_word(logp: float) -> str:
    if logp < 1:
        return "hydrophilic"
    if logp < 3:
        return "moderately lipophilic"
    if logp <= 5:
        return "lipophilic"
    return "highly lipophilic"


def summarize(druglikeness: Dict[str, Any]) -> str:
    """Render a plain-language summary from a ``calculate_druglikeness`` dict.

    Descriptive only: it characterizes the molecule's profile, it does not
    predict biological behavior.
    """
    mw = druglikeness["molecular_weight"]
    logp = druglikeness["logp"]
    parts = [f"A {_size_word(mw)}, {_lipophilicity_word(logp)} molecule (MW {mw}, logP {logp})."]

    lip = druglikeness["lipinski"]["violations"]
    if not lip:
        parts.append("Satisfies Lipinski's rule of five with no violations.")
    elif len(lip) == 1:
        parts.append(f"Within Lipinski's rule of five (one violation: {lip[0]}).")
    else:
        parts.append(f"Outside Lipinski's rule of five ({len(lip)} violations: {', '.join(lip)}).")

    veber = druglikeness["veber"]["violations"]
    if not veber:
        parts.append(
            f"Meets the Veber criteria (TPSA {druglikeness['tpsa']}, "
            f"{druglikeness['rotatable_bonds']} rotatable bonds)."
        )
    else:
        parts.append(f"Outside the Veber criteria ({', '.join(veber)}).")

    qed = druglikeness["qed"]
    qed_word = "high" if qed >= 0.67 else "moderate" if qed >= 0.5 else "low"
    parts.append(f"QED {qed} ({qed_word} drug-likeness).")

    return " ".join(parts)


def interpret(smiles: str) -> Dict[str, Any]:
    """Convenience: ``calculate_druglikeness`` plus a ``summary`` string."""
    druglikeness = calculate_druglikeness(smiles)
    return {**druglikeness, "summary": summarize(druglikeness)}
