"""Unified RNA-rSM-Bench v0 assembler.

Concatenates all per-source parquets into a single train/val/test
benchmark, then materializes the four OOD splits documented in
:mod:`boltzrna_diff.data.splits`.

Sources merged (in order, missing-tolerant):

* ``data/processed/hariboss.parquet``               (StructureRecord)
* ``data/processed/robin_binary.parquet``           (BinaryBinderRecord)
* ``data/processed/robin_per_target.parquet``       (InteractionRecord)
* ``data/processed/smrtnet.parquet``                (InteractionRecord)
* ``data/processed/rnamigos2_binary.parquet``       (InteractionRecord)
* ``data/processed/rnamigos2_docking.parquet``      (InteractionRecord)

The merged frame keeps a ``source`` column verbatim, plus a
``record_kind`` column distinguishing structures vs interactions vs
binders. Splits are emitted **per record_kind** so that, e.g., a
scaffold split for the binder-classification task is never contaminated
by the docking-score regression task.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .splits import SplitResult, split as do_split

log = logging.getLogger(__name__)

DEFAULT_SOURCES: dict[str, tuple[str, str]] = {
    # name → (relative_path, record_kind)
    "hariboss": ("hariboss.parquet", "structure"),
    "robin_binary": ("robin_binary.parquet", "binder"),
    "robin_per_target": ("robin_per_target.parquet", "interaction"),
    "smrtnet": ("smrtnet.parquet", "interaction"),
    "rnamigos2_binary": ("rnamigos2_binary.parquet", "interaction"),
    "rnamigos2_docking": ("rnamigos2_docking.parquet", "interaction"),
}


@dataclass
class BenchmarkBuild:
    df: pd.DataFrame                       # the merged frame
    splits: dict[str, dict[str, SplitResult]]  # record_kind → method → SplitResult


def _safe_load(path: Path, kind: str) -> pd.DataFrame | None:
    if not path.exists():
        log.warning("missing %s — skipping", path)
        return None
    df = pd.read_parquet(path)
    df = df.copy()
    df["record_kind"] = kind
    log.info("loaded %s | rows=%d", path.name, len(df))
    return df


def merge_sources(
    processed_dir: Path,
    sources: dict[str, tuple[str, str]] | None = None,
) -> pd.DataFrame:
    """Concatenate per-source parquets into one frame with ``record_kind``."""

    sources = sources or DEFAULT_SOURCES
    parts: list[pd.DataFrame] = []
    for name, (rel, kind) in sources.items():
        df = _safe_load(processed_dir / rel, kind)
        if df is None:
            continue
        df["bench_source"] = name
        parts.append(df)
    if not parts:
        raise RuntimeError(f"no parquets found under {processed_dir}")
    merged = pd.concat(parts, ignore_index=True, sort=False)
    log.info("merged: %d rows | record_kinds=%s", len(merged), dict(merged["record_kind"].value_counts()))
    return merged


def build_all_splits(
    df: pd.DataFrame,
    *,
    methods: list[str] | None = None,
    fractions: tuple[float, float, float] = (0.8, 0.1, 0.1),
) -> dict[str, dict[str, SplitResult]]:
    """Materialize the requested OOD splits per ``record_kind``.

    Splits are computed independently per ``record_kind`` so different
    tasks (structure / binder / interaction) never contaminate one
    another's train-vs-test partitions.
    """

    methods = methods or ["scaffold", "temporal", "family", "binding_site"]
    out: dict[str, dict[str, SplitResult]] = {}
    for kind, sub in df.groupby("record_kind"):
        if "ligand_smiles" not in sub.columns or sub["ligand_smiles"].isna().all():
            log.info("[%s] skipping ligand-based splits (no SMILES column)", kind)
            continue
        kind_out: dict[str, SplitResult] = {}
        for m in methods:
            try:
                if m == "temporal":
                    kind_out[m] = do_split(
                        sub,
                        method=m,
                        val_fraction=fractions[1],
                        test_fraction=fractions[2],
                    )
                else:
                    kind_out[m] = do_split(sub, method=m, fractions=fractions)
                log.info("[%s/%s] sizes=%s", kind, m, kind_out[m].sizes)
            except Exception as exc:  # noqa: BLE001
                log.warning("[%s/%s] split failed: %s", kind, m, exc)
        out[kind] = kind_out
    return out


def assignments_to_dataframe(
    df: pd.DataFrame,
    splits: dict[str, dict[str, SplitResult]],
) -> pd.DataFrame:
    """Add per-method split columns (``split_<method>``) to ``df``."""

    df = df.copy()
    for method in {m for d in splits.values() for m in d}:
        col = f"split_{method}"
        df[col] = pd.NA
        for kind, methods in splits.items():
            if method not in methods:
                continue
            assigns = methods[method].assignments()
            mask = df["record_kind"] == kind
            df.loc[mask, col] = df.loc[mask, "record_id"].map(assigns)
    return df


def build(
    processed_dir: Path,
    out_path: Path,
    *,
    methods: list[str] | None = None,
    fractions: tuple[float, float, float] = (0.8, 0.1, 0.1),
) -> BenchmarkBuild:
    """End-to-end driver: merge, split, write parquet."""

    merged = merge_sources(processed_dir)
    splits = build_all_splits(merged, methods=methods, fractions=fractions)
    out_df = assignments_to_dataframe(merged, splits)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out_path)
    log.info("wrote %s | rows=%d", out_path, len(out_df))
    return BenchmarkBuild(df=out_df, splits=splits)
