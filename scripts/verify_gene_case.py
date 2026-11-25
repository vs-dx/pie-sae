"""Read-only validation of PIE gene-symbol case normalization.

Run from the repository root:
    python scripts/verify_gene_case.py

Only reads saved gene vocabulary and biological references; no model inference.
This is a preliminary structural check, not a replacement for rerunning 06a/06b.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VOCAB = ROOT / "results/gene_reference/gene_index_to_symbol.json"
REFERENCES = {
    "GO_BP": ROOT / "data/biological_references/go_bp_gene_sets.json",
    "KEGG": ROOT / "data/biological_references/kegg_gene_sets.json",
    "Reactome": ROOT / "data/biological_references/reactome_gene_sets.json",
}


def main() -> None:
    genes = json.loads(VOCAB.read_text(encoding="utf-8"))["genes"]
    assert len(genes) == 6642
    canonical = [str(g).upper() for g in genes]
    assert len(set(canonical)) == 6642, "Case-normalized symbols are not unique"
    print("PIE gene vocabulary:", len(genes))
    print("Symbols requiring case normalization:", sum(a != b for a, b in zip(genes, canonical)))
    universe = set(canonical)
    counts = {}
    eligible = {}
    for source, path in REFERENCES.items():
        raw = json.loads(path.read_text(encoding="utf-8"))
        assert isinstance(raw, dict), (source, type(raw))
        count = 0
        for term, names in raw.items():
            members = universe.intersection(str(name).upper() for name in names)
            if 5 <= len(members) <= 500:
                eligible[(source, str(term))] = members
                count += 1
        counts[source] = count
    print("Eligible ontology terms by source:", counts)
    print("Total case-normalized ontology terms:", len(eligible))
    assert len(eligible) == 4421, "Expected corrected reference family of 4421 ontology terms"

    examples = [(key, members) for key, members in eligible.items()
                if "GO:0010648" in key[1]]
    assert len(examples) == 1, "Expected the recovered GO:0010648 term"
    assert "C6ORF89" in examples[0][1], "C6orf89 did not match C6ORF89"
    print("GO:0010648 boundary case: PASS (C6ORF89 present)")
    print("CASE NORMALIZATION PRECHECK: PASS")
    print("No original data or model artifacts were modified.")


if __name__ == "__main__":
    main()
