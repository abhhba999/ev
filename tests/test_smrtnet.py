"""Tests for the SMRTnet loader."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from boltzrna_diff.data.schema import InteractionRecord
from boltzrna_diff.data.smrtnet import _BENCHMARK_FILES, _TRAINING_FILE, load


RAW = Path("data/raw/smrtnet")


def _present() -> bool:
    return RAW.exists() and (RAW / _TRAINING_FILE).exists()


@pytest.mark.skipif(not _present(), reason="raw SMRTnet TSVs not downloaded")
def test_smrtnet_train_shape():
    out = load(RAW)
    df = out["smrtnet_train"]
    assert len(df) == 24687, "primary training set should be 24,687 rows"
    assert df["is_active"].sum() > 0
    assert (df["rna_length"] == 31).mean() > 0.99, "training rows are 31-nt sliding windows"


@pytest.mark.skipif(not _present(), reason="raw SMRTnet TSVs not downloaded")
def test_smrtnet_benchmarks_loaded():
    out = load(RAW)
    for tag in _BENCHMARK_FILES:
        assert tag in out, f"missing benchmark: {tag}"
        assert len(out[tag]) > 0, f"empty benchmark: {tag}"


@pytest.mark.skipif(not _present(), reason="raw SMRTnet TSVs not downloaded")
def test_smrtnet_record_validation():
    out = load(RAW)
    df = out["smrtnet_train"].head(200)
    valid = 0
    for _, row in df.iterrows():
        try:
            InteractionRecord(**{k: v for k, v in row.items() if pd.notna(v)})
            valid += 1
        except Exception:  # noqa: BLE001
            pass
    assert valid >= 195, f"expected ≥195/200 valid, got {valid}"
