"""Smoke tests for the ROBIN loader (read processed parquet, not raw)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

PARQUET = Path("data/processed/robin_per_target.parquet")
BIN_PARQUET = Path("data/processed/robin_binary.parquet")


@pytest.mark.skipif(not PARQUET.exists(), reason="run `boltzrna-diff data build-robin` first")
def test_robin_per_target_shape():
    df = pd.read_parquet(PARQUET)
    assert len(df) > 100_000
    # 25 RNA targets in ROBIN
    assert df["rna_sequence"].nunique() == 25
    assert df["is_active"].sum() > 1000


@pytest.mark.skipif(not PARQUET.exists(), reason="run `boltzrna-diff data build-robin` first")
def test_robin_pre_miR_21_present():
    df = pd.read_parquet(PARQUET)
    pre_miR_21_seq = "GGGUUGACUGUUGAAUCUCAUGGCAACCC"
    rows = df[df["rna_sequence"] == pre_miR_21_seq]
    assert len(rows) > 1000  # ROBIN screened ~24k compounds against pre-miR-21
    assert int(rows["is_active"].sum()) > 0


@pytest.mark.skipif(not BIN_PARQUET.exists(), reason="run `boltzrna-diff data build-robin` first")
def test_robin_binary_record_validation():
    from boltzrna_diff.data.robin import to_binary_records

    df = pd.read_parquet(BIN_PARQUET)
    # Validate first 200 rows through Pydantic.
    n = 0
    for r in to_binary_records(df.head(200)):
        n += 1
        assert r.is_rna_binder in (True, False)
    assert n >= 195  # >97% pass
