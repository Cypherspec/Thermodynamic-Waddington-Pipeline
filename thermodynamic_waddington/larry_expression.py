from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from typing import Any, Sequence


def _read_metadata(path: Path, max_cells: int | None = None) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader) if max_cells is None else [row for _, row in zip(range(max_cells), reader)]


def _read_gene_names(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


def _read_matrix_market_header(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError("not a Matrix Market file")
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        return tuple(map(int, line.split()))


def build_larry_expression_bundle(root: str | Path = "data/real/larry", max_cells: int | None = None, gene_indices: Sequence[int] | None = None, timepoint: float | None = 2.0) -> dict[str, Any]:
    """Build an aligned sparse expression bundle without densifying 130k cells."""
    root = Path(root)
    metadata = _read_metadata(root / "stateFate_inVitro_metadata.txt.gz", max_cells)
    if timepoint is not None:
        metadata = [row for row in metadata if float(row["Time point"]) == float(timepoint)]
    rows, genes, nonzeros = _read_matrix_market_header(root / "stateFate_inVitro_normed_counts.mtx.gz")
    gene_names = _read_gene_names(root / "stateFate_inVitro_gene_names.txt.gz")
    selected_genes = list(gene_indices) if gene_indices is not None else list(range(min(256, genes)))
    selected_set = set(selected_genes)
    cell_ids = [f"{row['Library']}::{row['Cell barcode']}" for row in metadata]
    labels = [row["Cell type annotation"] for row in metadata]
    coordinates = [[float(row["SPRING-x"]), float(row["SPRING-y"])] for row in metadata]
    return {
        "expression": [],
        "velocity": [],
        "labels": labels,
        "cell_ids": cell_ids,
        "embedding": coordinates,
        "gene_names": [gene_names[index] for index in selected_genes if index < len(gene_names)],
        "metadata": {
            "dataset": "Weinreb2020_LARRY_in_vitro",
            "timepoint": timepoint,
            "rows_in_source": rows,
            "genes_in_source": genes,
            "nonzero_entries": nonzeros,
            "selected_genes": selected_genes,
            "expression_storage": "sparse_matrix_market",
            "velocity_status": "not_observed",
            "claim_boundary": "LARRY metadata and coordinates are aligned; expression rows require sparse streaming before fitting.",
        },
        "source_schema": {"cells": rows, "genes": genes, "nonzeros": nonzeros, "selected_gene_count": len(selected_set)},
    }


def stream_larry_matrix(root: str | Path = "data/real/larry", selected_cells: Sequence[int] | None = None, selected_genes: Sequence[int] | None = None):
    """Yield sparse Matrix Market triplets as (cell, gene, value)."""
    root = Path(root)
    cell_filter = set(selected_cells) if selected_cells is not None else None
    gene_filter = set(selected_genes) if selected_genes is not None else None
    with gzip.open(root / "stateFate_inVitro_normed_counts.mtx.gz", "rt") as handle:
        handle.readline()
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            gene, cell, value = line.split()
            gene_index = int(gene) - 1
            cell_index = int(cell) - 1
            if (cell_filter is None or cell_index in cell_filter) and (gene_filter is None or gene_index in gene_filter):
                yield cell_index, gene_index, float(value)


def write_manifest(path: str | Path = "experiments/larry_expression_manifest.json") -> dict[str, Any]:
    root = Path("data/real/larry")
    bundle = build_larry_expression_bundle(root, max_cells=1)
    manifest = {"dataset": bundle["metadata"]["dataset"], "source": "https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/", "schema": bundle["source_schema"], "timepoints": [2.0, 4.0, 6.0], "storage": bundle["metadata"]["expression_storage"], "streaming_api": "thermodynamic_waddington.larry_expression.stream_larry_matrix", "status": "metadata_and_sparse_source_verified"}
    Path(path).write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
