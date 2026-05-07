"""OOD splits for RNA-rSM-Bench.

Four split strategies are supported, each addressing a distinct generalization
question:

1. **scaffold**  — split by Bemis-Murcko scaffold of the ligand. Tests
   whether the model generalizes beyond its training scaffolds. This is the
   default split used in PDBbind-style benchmarks.
2. **temporal**  — split by deposit / publication date. Holds out the most
   recent fraction. Mirrors a real prospective scenario (train on the past,
   evaluate on the future). Most relevant for HARIBOSS / R-BIND.
3. **family**    — split by RNA structural family / `rna_class`. Tests
   transfer from one RNA family to another (e.g. learn on riboswitches,
   evaluate on pre-miRNAs).
4. **binding_site** — split by RNA pocket / sequence cluster (CD-HIT style at
   ≥80 % identity over the binding-site residues). Tests whether the model
   has learned chemistry-specific binding rather than memorizing pocket
   sequences. Falls back to whole-sequence clustering when explicit
   pocket residue annotations are unavailable.

All splits return three disjoint sets of `record_id`s. Sizes default to
80/10/10 train/val/test but are tunable.
"""

from __future__ import annotations

import hashlib
import logging
import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

log = logging.getLogger(__name__)

SplitName = Literal["train", "val", "test"]


@dataclass
class SplitResult:
    """Container for a deterministic 3-way split."""

    train: list[str]
    val: list[str]
    test: list[str]

    def as_dict(self) -> dict[str, list[str]]:
        return {"train": self.train, "val": self.val, "test": self.test}

    @property
    def sizes(self) -> dict[str, int]:
        return {k: len(v) for k, v in self.as_dict().items()}

    def assignments(self) -> dict[str, SplitName]:
        a: dict[str, SplitName] = {}
        for split, ids in self.as_dict().items():
            for i in ids:
                a[i] = split  # type: ignore[assignment]
        return a


# ----- helpers -----

def _scaffold_smiles(smi: str) -> str:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return ""
    try:
        scaff = MurckoScaffold.GetScaffoldForMol(m)
    except Exception:  # noqa: BLE001
        return ""
    return Chem.MolToSmiles(scaff, canonical=True)


def _stable_bucket(key: str, total: int) -> int:
    """Hash a key to an integer bucket in ``[0, total)`` deterministically."""

    h = hashlib.blake2b(key.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(h, byteorder="big") % total


def _greedy_pack(
    grouped_ids: dict[str, list[str]],
    fractions: tuple[float, float, float],
    seed: int,
) -> SplitResult:
    """Pack groups into 3 buckets by descending size (deterministic)."""

    rng = np.random.default_rng(seed)
    keys = sorted(grouped_ids, key=lambda k: (-len(grouped_ids[k]), k))
    n_total = sum(len(v) for v in grouped_ids.values())
    targets = [int(round(f * n_total)) for f in fractions]
    targets[-1] = n_total - sum(targets[:-1])  # absorb rounding into test
    bucket_ids: list[list[str]] = [[], [], []]
    bucket_sz = [0, 0, 0]
    for k in keys:
        # pick least-full bucket relative to its target capacity
        order = sorted(range(3), key=lambda i: bucket_sz[i] / max(targets[i], 1))
        # tie-break with deterministic random
        order_arr = np.array(order)
        rng.shuffle(order_arr)
        chosen = int(order_arr[0])
        bucket_ids[chosen].extend(grouped_ids[k])
        bucket_sz[chosen] += len(grouped_ids[k])
    return SplitResult(train=bucket_ids[0], val=bucket_ids[1], test=bucket_ids[2])


# ----- 1. scaffold split -----

def scaffold_split(
    df: pd.DataFrame,
    *,
    smiles_col: str = "ligand_smiles",
    id_col: str = "record_id",
    fractions: tuple[float, float, float] = (0.8, 0.1, 0.1),
    seed: int = 0,
) -> SplitResult:
    """Bemis-Murcko scaffold split. Compounds without a valid scaffold (mostly
    chains) are bucketed by canonical SMILES instead.
    """

    grouped: dict[str, list[str]] = defaultdict(list)
    for _, r in df.iterrows():
        rid = str(r[id_col])
        scaff = _scaffold_smiles(str(r[smiles_col]))
        if not scaff:
            scaff = f"FALLBACK::{r[smiles_col]}"
        grouped[scaff].append(rid)
    result = _greedy_pack(grouped, fractions, seed)
    log.info("scaffold_split: %s scaffolds, sizes=%s", len(grouped), result.sizes)
    return result


# ----- 2. temporal split -----

def temporal_split(
    df: pd.DataFrame,
    *,
    date_col: str = "deposit_date",
    id_col: str = "record_id",
    val_fraction: float = 0.1,
    test_fraction: float = 0.1,
) -> SplitResult:
    """Sort by date (ISO strings sort lexicographically); newest val_fraction
    + test_fraction become val/test.
    """

    sub = df[[id_col, date_col]].dropna(subset=[date_col]).copy()
    sub[date_col] = sub[date_col].astype(str)
    sub = sub.sort_values(date_col)
    n = len(sub)
    n_test = int(round(n * test_fraction))
    n_val = int(round(n * val_fraction))
    n_train = n - n_test - n_val
    train = sub.iloc[:n_train][id_col].tolist()
    val = sub.iloc[n_train : n_train + n_val][id_col].tolist()
    test = sub.iloc[n_train + n_val :][id_col].tolist()
    # Records without dates: append to train (conservative).
    missing = df[~df[id_col].isin(sub[id_col])][id_col].tolist()
    train.extend(missing)
    log.info("temporal_split: n=%d, train=%d val=%d test=%d", n, len(train), len(val), len(test))
    return SplitResult(train=train, val=val, test=test)


# ----- 3. family split -----

def family_split(
    df: pd.DataFrame,
    *,
    family_col: str = "rna_class",
    id_col: str = "record_id",
    fractions: tuple[float, float, float] = (0.8, 0.1, 0.1),
    seed: int = 0,
    val_families: Iterable[str] | None = None,
    test_families: Iterable[str] | None = None,
) -> SplitResult:
    """Split by RNA family. Either pass explicit ``val_families`` /
    ``test_families`` to fix the holdout, or let the function pack families
    greedily into the three buckets by size.
    """

    grouped: dict[str, list[str]] = defaultdict(list)
    for _, r in df.iterrows():
        grouped[str(r[family_col])].append(str(r[id_col]))

    if val_families is not None or test_families is not None:
        train, val, test = [], [], []
        val_set = set(val_families or [])
        test_set = set(test_families or [])
        for fam, ids in grouped.items():
            if fam in test_set:
                test.extend(ids)
            elif fam in val_set:
                val.extend(ids)
            else:
                train.extend(ids)
        log.info(
            "family_split (explicit holdout): val_families=%s test_families=%s sizes=%s",
            val_set, test_set, (len(train), len(val), len(test)),
        )
        return SplitResult(train=train, val=val, test=test)

    result = _greedy_pack(grouped, fractions, seed)
    log.info("family_split (greedy): %d families sizes=%s", len(grouped), result.sizes)
    return result


# ----- 4. binding-site split -----

def binding_site_split(
    df: pd.DataFrame,
    *,
    seq_col: str = "rna_sequence",
    id_col: str = "record_id",
    pocket_col: str | None = "pocket_residues",
    similarity_threshold: float = 0.8,
    fractions: tuple[float, float, float] = (0.8, 0.1, 0.1),
    seed: int = 0,
) -> SplitResult:
    """Cluster by RNA sequence (or pocket-residue subsequence if available)
    using a length-bucketed shingle hash; compounds whose RNA targets
    cluster together end up in the same split.

    This is a *lightweight* approximation of CD-HIT clustering (no third-
    party install required for CPU CI). Two sequences cluster together when:

    - they have the same length, AND
    - they share ≥ ``similarity_threshold`` fraction of overlapping
      length-3 shingles.

    For HARIBOSS the resolution is sufficient (most binding sites are <50 nt
    and well-conserved within an RNA family); we will swap in MMseqs2 / CD-
    HIT-EST on the GPU cluster for production runs.
    """

    rng = np.random.default_rng(seed)

    def _shingles(s: str, k: int = 3) -> set[str]:
        return {s[i : i + k] for i in range(len(s) - k + 1)} if len(s) >= k else {s}

    if pocket_col and pocket_col in df.columns:
        cluster_seq = df[pocket_col].fillna(df[seq_col]).astype(str)
    else:
        cluster_seq = df[seq_col].astype(str)

    # Greedy single-link clustering over unique sequences (small support).
    unique_seqs = sorted(set(cluster_seq.tolist()))
    log.info("binding_site_split: %d unique sequences before clustering", len(unique_seqs))
    parents: dict[str, str] = {s: s for s in unique_seqs}

    def _find(s: str) -> str:
        while parents[s] != s:
            parents[s] = parents[parents[s]]
            s = parents[s]
        return s

    def _union(a: str, b: str) -> None:
        parents[_find(a)] = _find(b)

    # Index sequences by length to limit pairwise comparisons.
    by_length: dict[int, list[str]] = defaultdict(list)
    for s in unique_seqs:
        by_length[len(s)].append(s)
    for length, group in by_length.items():
        if length < 6:
            for s in group:
                _union(s, group[0])
            continue
        for i, a in enumerate(group):
            sa = _shingles(a)
            for b in group[i + 1 :]:
                sb = _shingles(b)
                jacc = len(sa & sb) / max(1, len(sa | sb))
                if jacc >= similarity_threshold:
                    _union(a, b)

    cluster_of = {s: _find(s) for s in unique_seqs}
    grouped: dict[str, list[str]] = defaultdict(list)
    for _, r in df.iterrows():
        rid = str(r[id_col])
        s = str(r.get(pocket_col, "") or r[seq_col]) if pocket_col else str(r[seq_col])
        c = cluster_of.get(s, s)
        grouped[c].append(rid)

    log.info("binding_site_split: %d clusters", len(grouped))
    return _greedy_pack(grouped, fractions, seed)


# ----- driver -----

_SPLITTERS = {
    "scaffold": scaffold_split,
    "temporal": temporal_split,
    "family": family_split,
    "binding_site": binding_site_split,
}


def split(
    df: pd.DataFrame,
    *,
    method: Literal["scaffold", "temporal", "family", "binding_site"],
    **kwargs,
) -> SplitResult:
    """Driver for split methods. ``method`` selects the splitter."""

    if method not in _SPLITTERS:
        raise ValueError(f"unknown split method '{method}', choose from {list(_SPLITTERS)}")
    return _SPLITTERS[method](df, **kwargs)
