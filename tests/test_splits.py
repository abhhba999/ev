"""Tests for OOD split strategies."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from boltzrna_diff.data.splits import (
    SplitResult,
    binding_site_split,
    family_split,
    scaffold_split,
    split as split_driver,
    temporal_split,
)


@pytest.fixture
def toy_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"record_id": "a", "ligand_smiles": "c1ccccc1", "rna_sequence": "GCUAGCUA",
             "rna_class": "riboswitch", "deposit_date": "2010-01-01"},
            {"record_id": "b", "ligand_smiles": "c1ccccc1C",  "rna_sequence": "GCUAGCUA",
             "rna_class": "riboswitch", "deposit_date": "2012-05-12"},
            {"record_id": "c", "ligand_smiles": "CCO", "rna_sequence": "AAACCCUU",
             "rna_class": "pre_mirna",  "deposit_date": "2018-09-01"},
            {"record_id": "d", "ligand_smiles": "OCCN", "rna_sequence": "AAACCCUU",
             "rna_class": "pre_mirna",  "deposit_date": "2019-03-15"},
            {"record_id": "e", "ligand_smiles": "c1ccncc1", "rna_sequence": "UAGCUAGCU",
             "rna_class": "viral_rna",  "deposit_date": "2024-04-01"},
            {"record_id": "f", "ligand_smiles": "c1ccncc1Cl", "rna_sequence": "UAGCUAGCU",
             "rna_class": "viral_rna",  "deposit_date": "2025-06-30"},
        ]
    )


def _disjoint(s: SplitResult) -> bool:
    sets = [set(s.train), set(s.val), set(s.test)]
    return all(len(a & b) == 0 for i, a in enumerate(sets) for b in sets[i + 1 :])


def _covers_all(df: pd.DataFrame, s: SplitResult) -> bool:
    return set(s.train + s.val + s.test) == set(df["record_id"])


def test_scaffold_split_disjoint_and_covers(toy_df):
    s = scaffold_split(toy_df, fractions=(0.5, 0.25, 0.25), seed=0)
    assert _disjoint(s)
    assert _covers_all(toy_df, s)


def test_temporal_split_orders_by_date(toy_df):
    s = temporal_split(toy_df, val_fraction=0.0, test_fraction=1 / 3)
    # Newest 2 records should be in test.
    assert set(s.test) == {"e", "f"}
    assert _disjoint(s)


def test_family_split_explicit_holdout(toy_df):
    s = family_split(toy_df, val_families=["pre_mirna"], test_families=["viral_rna"])
    assert set(s.train) == {"a", "b"}
    assert set(s.val) == {"c", "d"}
    assert set(s.test) == {"e", "f"}


def test_binding_site_split_clusters_identical_sequences(toy_df):
    s = binding_site_split(toy_df, seq_col="rna_sequence", pocket_col=None,
                            fractions=(0.5, 0.25, 0.25), seed=0)
    # All compounds against same RNA must end up in the same bucket.
    assignments = s.assignments()
    assert assignments["a"] == assignments["b"]
    assert assignments["c"] == assignments["d"]
    assert assignments["e"] == assignments["f"]
    assert _disjoint(s)
    assert _covers_all(toy_df, s)


def test_split_driver_dispatch(toy_df):
    s = split_driver(toy_df, method="scaffold", fractions=(0.5, 0.25, 0.25), seed=0)
    assert _disjoint(s)


# ---- end-to-end on real ROBIN per-target parquet (skip if not built) ----

ROBIN_PARQUET = Path("data/processed/robin_per_target.parquet")


@pytest.mark.skipif(not ROBIN_PARQUET.exists(), reason="run `boltzrna-diff data build-robin` first")
def test_family_split_on_robin():
    df = pd.read_parquet(ROBIN_PARQUET)
    s = family_split(
        df,
        val_families=["pre_mirna"],
        test_families=["g_quadruplex"],
    )
    assert set(s.train) & set(s.val) == set()
    assert set(s.train) & set(s.test) == set()
    assert len(s.train) + len(s.val) + len(s.test) == len(df)


@pytest.mark.skipif(not ROBIN_PARQUET.exists(), reason="run `boltzrna-diff data build-robin` first")
def test_binding_site_split_on_robin():
    df = pd.read_parquet(ROBIN_PARQUET).head(20_000)  # subsample for speed
    s = binding_site_split(df, seq_col="rna_sequence", pocket_col=None,
                            fractions=(0.7, 0.15, 0.15), seed=0)
    # ROBIN's 25 RNA targets are all very different, so we expect ~25 clusters.
    sizes = s.sizes
    assert sum(sizes.values()) == len(df)
    # Each train/val/test bucket should have non-trivial size.
    assert all(v > 0 for v in sizes.values())
