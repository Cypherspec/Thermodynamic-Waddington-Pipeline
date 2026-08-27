from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any


PAUL15_URL = "https://exampledata.scverse.org/scanpy/paul15.h5"
PANCREAS_URL = "https://github.com/theislab/scvelo_notebooks/raw/master/data/Pancreas/endocrinogenesis_day15.h5ad"
LARRY_FIGSHARE_DOI = "https://doi.org/10.6084/m9.figshare.28491017.v1"
LARRY_FIGSHARE_PAGE = "https://figshare.com/articles/dataset/adata_Weinreb2020_in_vitro_gene_filtered_h5ad/28491017"


@dataclass(frozen=True)
class DatasetManifest:
    name: str
    source_url: str
    local_path: str
    expected_format: str
    download_status: str
    notes: str

    def to_dict(self) -> dict[str, str]:
        return asdict(self)


def dataset_catalog() -> list[DatasetManifest]:
    return [
        DatasetManifest("paul15", PAUL15_URL, "data/real/paul15.h5", "scanpy HDF5", "available-on-demand", "Canonical 2,730-cell Paul et al. hematopoiesis benchmark; raw count matrix and cluster annotations."),
        DatasetManifest("pancreas", PANCREAS_URL, "data/real/endocrinogenesis_day15.h5ad", "AnnData H5AD", "available-on-demand", "scVelo pancreas velocity benchmark; use only with spliced/unspliced layers or a separately generated velocity field."),
        DatasetManifest("larry_in_vitro", LARRY_FIGSHARE_DOI, "data/real/larry_in_vitro.h5ad", "AnnData H5AD", "metadata-only", "Lineage-traced Weinreb 2020 derivative. Select the concrete asset from the Figshare page before downloading; the project never guesses a large file URL."),
    ]


def catalog_payload() -> dict[str, object]:
    return {"datasets": [item.to_dict() for item in dataset_catalog()], "scientific_policy": "Real-data downloads are opt-in; no dataset is silently fetched during fitting.", "larry_figshare_page": LARRY_FIGSHARE_PAGE}


def write_catalog(path: str | Path) -> None:
    Path(path).write_text(json.dumps(catalog_payload(), indent=2, sort_keys=True))


def download_dataset(name: str, destination: str | Path, timeout: int = 120) -> DatasetManifest:
    matches = [item for item in dataset_catalog() if item.name == name]
    if not matches:
        raise ValueError(f"unknown dataset {name!r}; choose from {[item.name for item in dataset_catalog()]}")
    item = matches[0]
    if name == "larry_in_vitro":
        raise ValueError("LARRY requires selecting a concrete Figshare asset; use the DOI in the manifest rather than guessing a large file URL")
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(item.source_url, headers={"User-Agent": "thermodynamic-waddington/0.2"})
    with urllib.request.urlopen(request, timeout=timeout) as response, destination.open("wb") as handle:
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            handle.write(block)
    return DatasetManifest(item.name, item.source_url, str(destination), item.expected_format, "downloaded", item.notes)


def _try_import(module: str) -> Any:
    try:
        return __import__(module)
    except ImportError as error:
        raise RuntimeError(f"{module} is required for this workflow; install the optional real-data dependencies") from error


def inspect_paul15(path: str | Path) -> dict[str, object]:
    import h5py
    with h5py.File(path, "r") as handle:
        matrix = handle["data.debatched"]
        return {"path": str(path), "format": "Paul15 HDF5", "cells": int(matrix.shape[1]), "genes": int(matrix.shape[0]), "obs_columns": ["cluster.id"], "layers": [], "obsm": [], "velocity_status": "not_observed", "gene_name_dataset": "data.debatched_rownames"}


def inspect_anndata(path: str | Path, max_cells: int = 256, max_genes: int = 128) -> dict[str, object]:
    if str(path).lower().endswith((".h5", ".hdf5")):
        return inspect_paul15(path)
    anndata = _try_import("anndata")
    adata = anndata.read_h5ad(path, backed="r")
    obs_names = [str(item) for item in adata.obs_names[: min(max_cells, adata.n_obs)]]
    var_names = [str(item) for item in adata.var_names[: min(max_genes, adata.n_vars)]]
    return {"path": str(path), "cells": int(adata.n_obs), "genes": int(adata.n_vars), "obs_columns": [str(item) for item in adata.obs.columns], "var_columns": [str(item) for item in adata.var.columns], "layers": [str(item) for item in adata.layers.keys()], "obsm": [str(item) for item in adata.obsm.keys()], "obs_preview": obs_names[:8], "gene_preview": var_names[:8]}


def derive_velocity_from_spliced_unspliced(path: str | Path, max_cells: int | None = None, max_genes: int | None = None, label_key: str | None = None, layer_spliced: str = "spliced", layer_unspliced: str = "unspliced") -> dict[str, object]:
    anndata = _try_import("anndata")
    adata = anndata.read_h5ad(path)
    if max_cells is not None:
        adata = adata[:max_cells].copy()
    if max_genes is not None:
        adata = adata[:, :max_genes].copy()
    if layer_spliced not in adata.layers or layer_unspliced not in adata.layers:
        raise ValueError(f"both {layer_spliced!r} and {layer_unspliced!r} layers are required; available={list(adata.layers.keys())}")
    spliced = _dense(adata.layers[layer_spliced])
    unspliced = _dense(adata.layers[layer_unspliced])
    velocity = [[unspliced[i][j] - spliced[i][j] for j in range(len(spliced[i]))] for i in range(len(spliced))]
    labels = [str(value) for value in adata.obs[label_key]] if label_key and label_key in adata.obs else ["unlabeled"] * len(spliced)
    return {"expression": spliced, "velocity": velocity, "labels": labels, "cell_ids": [str(item) for item in adata.obs_names], "gene_names": [str(item) for item in adata.var_names], "metadata": {"source": str(path), "expression_key": layer_spliced, "velocity_key": "derived:unspliced-minus-spliced", "velocity_status": "derived_proxy", "derivation": "per-cell unspliced minus spliced; not scVelo dynamical-model velocity", "label_key": label_key}}


def _dense(value: Any) -> list[list[float]]:
    if hasattr(value, "toarray"):
        value = value.toarray()
    if hasattr(value, "tolist"):
        value = value.tolist()
    return [[float(item) for item in row] for row in value]


def anndata_to_bundle(path: str | Path, expression_key: str = "X", velocity_key: str = "velocity", label_key: str | None = None, max_cells: int | None = None, max_genes: int | None = None) -> dict[str, object]:
    if str(path).lower().endswith((".h5", ".hdf5")):
        return paul15_to_bundle(path, max_cells=max_cells, max_genes=max_genes)
    anndata = _try_import("anndata")
    adata = anndata.read_h5ad(path)
    if max_cells is not None:
        adata = adata[:max_cells].copy()
    if max_genes is not None:
        adata = adata[:, :max_genes].copy()
    expression = _dense(adata.X if expression_key == "X" else adata.layers[expression_key])
    if velocity_key in adata.layers:
        velocity = _dense(adata.layers[velocity_key])
    elif velocity_key in adata.obsm:
        velocity = _dense(adata.obsm[velocity_key])
    else:
        raise ValueError(f"velocity key {velocity_key!r} not found in layers or obsm; available layers={list(adata.layers.keys())}, obsm={list(adata.obsm.keys())}")
    labels = [str(value) for value in adata.obs[label_key]] if label_key and label_key in adata.obs else ["unlabeled"] * len(expression)
    return {"expression": expression, "velocity": velocity, "labels": labels, "cell_ids": [str(item) for item in adata.obs_names], "metadata": {"source": str(path), "expression_key": expression_key, "velocity_key": velocity_key, "label_key": label_key, "layers": list(adata.layers.keys()), "obsm": list(adata.obsm.keys()), "gene_names": [str(item) for item in adata.var_names]}}


def compute_file_sha256(path: str | Path, block_size: int = 1024 * 1024) -> str:
    import hashlib
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def benchmark_provenance(path: str | Path, dataset_name: str, source_url: str, observed_velocity: bool, expression_key: str = "X", velocity_key: str = "velocity", label_key: str | None = None) -> dict[str, object]:
    inspection = inspect_anndata(path)
    return {
        "dataset": dataset_name,
        "source_url": source_url,
        "local_path": str(path),
        "sha256": compute_file_sha256(path),
        "bytes": Path(path).stat().st_size,
        "schema": inspection,
        "representation": {"expression_key": expression_key, "velocity_key": velocity_key, "label_key": label_key},
        "velocity_observed": bool(observed_velocity),
        "claim_boundary": "Benchmark provenance and external-data reproducibility only; not biological validation by itself.",
    }


def validate_velocity_schema(path: str | Path, velocity_key: str = "velocity") -> dict[str, object]:
    if str(path).lower().endswith((".h5", ".hdf5")):
        return {"valid": False, "velocity_key": velocity_key, "reason": "Paul15 contains no measured RNA velocity"}
    anndata = _try_import("anndata")
    adata = anndata.read_h5ad(path, backed="r")
    candidates = list(adata.layers.keys()) + list(adata.obsm.keys())
    found = velocity_key in candidates
    return {"valid": found, "velocity_key": velocity_key, "available": candidates, "shape": [int(adata.n_obs), int(adata.n_vars)] if found and velocity_key in adata.layers else None, "reason": None if found else "Requested velocity representation is absent"}


def verify_download(path: str | Path, expected_format: str | None = None) -> dict[str, object]:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise ValueError(f"download is empty or missing: {path}")
    suffix = path.suffix.lower()
    recognized = suffix in {".h5", ".h5ad", ".loom", ".csv", ".tsv", ".json"}
    return {"path": str(path), "bytes": path.stat().st_size, "suffix": suffix, "recognized_suffix": recognized, "expected_format": expected_format}


def paul15_to_bundle(path: str | Path, max_cells: int | None = None, max_genes: int | None = None) -> dict[str, object]:
    try:
        import h5py
    except ImportError as error:
        raise RuntimeError("h5py is required to read Paul15 HDF5 data") from error
    with h5py.File(path, "r") as handle:
        matrix = handle["data.debatched"][()]
        genes = [str(value.decode() if hasattr(value, "decode") else value).split(";")[0] for value in handle["data.debatched_rownames"][()]]
        cell_ids = [str(value.decode() if hasattr(value, "decode") else value) for value in handle["data.debatched_colnames"][()]]
        clusters = [str(int(value[0] if hasattr(value, "__len__") else value)) for value in handle["cluster.id"][()]]
    selected_genes = list(range(min(len(genes), max_genes))) if max_genes is not None else list(range(len(genes)))
    selected_cells = list(range(min(len(cell_ids), max_cells))) if max_cells is not None else list(range(len(cell_ids)))
    expression = [[float(matrix[gene, cell]) for gene in selected_genes] for cell in selected_cells]
    labels = [clusters[cell] for cell in selected_cells]
    velocity = [[0.0 for _ in selected_genes] for _ in selected_cells]
    return {"expression": expression, "velocity": velocity, "labels": labels, "cell_ids": [cell_ids[cell] for cell in selected_cells], "gene_names": [genes[gene] for gene in selected_genes], "metadata": {"source": str(path), "dataset": "Paul15", "velocity_status": "not_observed", "warning": "Paul15 contains expression and clusters but no RNA velocity. Zero velocity is a placeholder and must not be interpreted as measured dynamics."}}
