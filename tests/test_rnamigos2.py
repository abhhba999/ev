"""Tests for the RNAmigos2 loader."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from boltzrna_diff.data.rnamigos2 import parse_pocket_json


RAW = Path("data/raw/rnamigos2")
POCKETS_DIR = RAW / "json_pockets_expanded"


def _pockets_present() -> bool:
    return POCKETS_DIR.exists() and any(POCKETS_DIR.glob("*.json"))


@pytest.fixture
def synthetic_pocket(tmp_path: Path) -> Path:
    """Tiny in-memory NetworkX pocket: 5-nt RNA, 2 in-pocket residues."""

    g = {
        "directed": True,
        "multigraph": False,
        "graph": {},
        "nodes": [
            {"in_pocket": False, "nt_code": "G", "id": "1xyz.A.1"},
            {"in_pocket": True,  "nt_code": "C", "id": "1xyz.A.2"},
            {"in_pocket": True,  "nt_code": "A", "id": "1xyz.A.3"},
            {"in_pocket": False, "nt_code": "U", "id": "1xyz.A.4"},
            {"in_pocket": False, "nt_code": "G", "id": "1xyz.A.5"},
        ],
        "links": [],
    }
    p = tmp_path / "1XYZ_A_LIG_1.json"
    p.write_text(json.dumps(g))
    return p


def test_parse_pocket_json_synthetic(synthetic_pocket):
    info = parse_pocket_json(synthetic_pocket)
    assert info is not None
    assert info.pdb_id == "1xyz"
    assert info.chain == "A"
    assert info.sequence == "GCAUG"
    assert info.pocket_residues == [1, 2]


def test_parse_pocket_json_handles_modified_residues(tmp_path: Path):
    g = {
        "directed": True,
        "multigraph": False,
        "graph": {},
        "nodes": [
            {"in_pocket": False, "nt_code": "G", "id": "1abc.A.1"},
            {"in_pocket": False, "nt_code": "X", "id": "1abc.A.2"},  # modified base, drop
            {"in_pocket": True, "nt_code": "C", "id": "1abc.A.3"},
        ],
        "links": [],
    }
    p = tmp_path / "1ABC_A_LIG_1.json"
    p.write_text(json.dumps(g))
    info = parse_pocket_json(p)
    assert info is not None
    assert info.sequence == "GC"
    assert info.pocket_residues == [1]


@pytest.mark.skipif(not _pockets_present(), reason="rnamigos2 pockets not extracted")
def test_real_rnamigos_pockets_load():
    from boltzrna_diff.data.rnamigos2 import load_pockets

    pockets = load_pockets(RAW)
    assert len(pockets) > 1500, "expected at least 1,500 RNAmigos2 pockets"
    seq_lens = [len(p.sequence) for p in pockets.values()]
    assert max(seq_lens) > 50
    assert min(seq_lens) >= 1
