"""Prepare STRING v12.0 human functional-association pairs for PIE gene symbols.

Source files are downloaded only once; compressed source files and their hashes
are kept. STRING identifiers are mapped through the corresponding v12.0
protein.info preferred_name column. The resulting edges are undirected.

Reproduction note: STRING's 'protein.links' is an association network, not
restricted to experimentally verified physical protein interactions.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
from urllib.request import Request, urlopen

VERSION = "12.0"
BASE_URL = "https://stringdb-downloads.org/download"
FILES = {
    "links": (
        "9606.protein.links.v12.0.txt.gz",
        "protein.links.v12.0",
    ),
    "info": (
        "9606.protein.info.v12.0.txt.gz",
        "protein.info.v12.0",
    ),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(2 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_if_absent(dest: Path, url: str) -> None:
    if dest.is_file():
        if dest.stat().st_size < 1024:
            raise ValueError(f"Existing source file is unexpectedly small: {dest}")
        print(f"Using existing download: {dest.name}")
        return

    temp = dest.with_name(dest.name + ".partial")
    if temp.exists():
        temp.unlink()  # previous incomplete download only

    print(f"Downloading {dest.name} ...", flush=True)
    request = Request(url, headers={"User-Agent": "PIE-SAE-research/1.0"})
    try:
        with urlopen(request, timeout=120) as response, temp.open("wb") as out:
            for chunk in iter(lambda: response.read(1024 * 1024), b""):
                out.write(chunk)
        with gzip.open(temp, "rb") as check:
            if not check.read(128):
                raise ValueError(f"Empty compressed source: {dest.name}")
        temp.replace(dest)
    except Exception:
        if temp.exists():
            temp.unlink()
        raise
    print(f"Downloaded {dest.name} ({dest.stat().st_size / (1024**2):.1f} MiB)")


def prepare_string_reference(
    project_root: str | Path,
    *,
    reference_dir: str | Path | None = None,
    source_paths: dict[str, str | Path] | None = None,
) -> dict:
    """Prepare a JSON containing STRING gene-symbol edges at scores >=700, >=900.

    source_paths can override {"info": path, "links": path} for offline testing.
    """
    project_root = Path(project_root).resolve()
    reference_dir = (
        Path(reference_dir)
        if reference_dir is not None
        else project_root / "data" / "biological_references"
    )
    reference_dir.mkdir(parents=True, exist_ok=True)

    vocab_path = project_root / "results" / "gene_reference" / "gene_index_to_symbol.json"
    with vocab_path.open(encoding="utf-8") as f:
        vocabulary = json.load(f)
    pie_genes = {str(g).upper() for g in vocabulary["genes"]}
    assert len(pie_genes) == 6642, "Expected 6642 unique PIE gene symbols"

    sources = {}
    for key, (filename, url_folder) in FILES.items():
        path = (
            Path(source_paths[key])
            if source_paths is not None
            else reference_dir / filename
        )
        url = f"{BASE_URL}/{url_folder}/{filename}"
        if source_paths is None:
            download_if_absent(path, url)
        elif not path.is_file():
            raise FileNotFoundError(path)
        sources[key] = {"path": path, "url": url, "sha256": sha256_file(path)}

    # Map STRING proteins to their human gene symbols; some symbols may identify
    # more than one protein. The final network is collapsed to distinct gene pairs.
    protein_to_symbol = {}
    with gzip.open(sources["info"]["path"], "rt", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        if len(header) < 2 or "preferred_name" not in header:
            raise ValueError(f"Unexpected protein.info header: {header}")
        name_idx = header.index("preferred_name")
        id_idx = 0
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) <= name_idx:
                continue
            name = fields[name_idx].strip().upper()
            if name in pie_genes:
                protein_to_symbol[fields[id_idx]] = name

    mapped_genes = set(protein_to_symbol.values())
    print(f"PIE genes with matching STRING preferred name: {len(mapped_genes)}/6642")
    if len(mapped_genes) < 100:  # deliberately liberal: catch schema mismatch
        raise ValueError("Very few mapped genes; check STRING species and info header")

    edges_700: set[tuple[str, str]] = set()
    edges_900: set[tuple[str, str]] = set()
    scanned = 0
    with gzip.open(sources["links"]["path"], "rt", encoding="utf-8") as f:
        header = f.readline().split()
        if header[:3] != ["protein1", "protein2", "combined_score"]:
            raise ValueError(f"Unexpected protein.links header: {header}")

        for line in f:
            scanned += 1
            left, right, score = line.split()[:3]
            score = int(score)
            if score < 700:
                continue
            a = protein_to_symbol.get(left)
            b = protein_to_symbol.get(right)
            if not a or not b or a == b:
                continue
            pair = (a, b) if a < b else (b, a)
            edges_700.add(pair)
            if score >= 900:
                edges_900.add(pair)

    assert edges_900.issubset(edges_700)
    if not edges_700:
        raise ValueError("No STRING associations remain after filtering")

    payload = {
        "pairs_700": [list(pair) for pair in sorted(edges_700)],
        "pairs_900": [list(pair) for pair in sorted(edges_900)],
    }
    output_path = reference_dir / "string_ppi_edges.json"
    if output_path.exists():
        with output_path.open(encoding="utf-8") as f:
            old = json.load(f)
        if old != payload:
            raise FileExistsError(
                f"Existing STRING network differs; refusing overwrite: {output_path}"
            )
        print("Existing output matches generated network; no overwrite.")
    else:
        tmp_path = output_path.with_suffix(".json.partial")
        with tmp_path.open("w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False)
        os.replace(tmp_path, output_path)

    result = {
        "database": "STRING",
        "version": VERSION,
        "species": "Homo sapiens (9606)",
        "network": "functional protein associations; not exclusively physical interactions",
        "identifier_mapping": "protein.info preferred_name -> PIE gene symbol (case-insensitive)",
        "thresholds": [700, 900],
        "n_gene_symbols_mapped": len(mapped_genes),
        "n_source_links_scanned": scanned,
        "n_pairs_700": len(edges_700),
        "n_pairs_900": len(edges_900),
        "output_file": output_path.name,
        "source_files": {
            key: {
                "url": item["url"],
                "file": item["path"].name,
                "sha256": item["sha256"],
            }
            for key, item in sources.items()
        },
    }

    provenance_path = reference_dir / "provenance.json"
    if provenance_path.is_file():
        with provenance_path.open(encoding="utf-8") as f:
            provenance = json.load(f)
    else:
        provenance = {"sources": {}}
    provenance.setdefault("sources", {})["STRING"] = result
    with provenance_path.open("w", encoding="utf-8") as f:
        json.dump(provenance, f, indent=2, ensure_ascii=False)

    print("STRING v12.0 human network prepared")
    print(f"Gene pairs at score >=700: {len(edges_700)}")
    print(f"Gene pairs at score >=900: {len(edges_900)}")
    print(f"Output: {output_path}")
    print(f"Provenance: {provenance_path}")
    return result
