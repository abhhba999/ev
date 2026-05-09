"""SMRTnet loader.

SMRTnet (Fei et al., *Nature Biotechnology* 2026, doi:10.1038/s41587-025-02942-z)
is a deep-learning method for predicting RNA-small-molecule interactions from
RNA secondary structure. The repository ships with a primary training set
(``SMRTnet_data.txt``) and five literature-derived benchmarks (R-BIND, NALDB,
SMMRNA, RSIM, NewPub).

GitHub: https://github.com/Yuhan-Fei/SMRTnet (data/ folder)

File formats
------------
- ``SMRTnet_data.txt`` — 24,687 rows × 4 cols::

      SMILES \t RNA_seq(31nt) \t dot_bracket(31) \t label

- ``SMRTnet_benchmark*.txt`` — 5 cols::

      id \t SMILES \t RNA_seq \t dot_bracket \t label

- ``MYC_IRES.txt`` / ``MYC_RIBOTAC.txt`` / ``natural_compounds.txt`` — small
  utility tables (skipped by this loader).

Each row maps to a :class:`~boltzrna_diff.data.schema.InteractionRecord`,
where the binary `label` becomes ``is_active``.
"""

from __future__ import annotations

import logging
import re
import subprocess
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from tqdm import tqdm

from .schema import AssayKind, InteractionRecord, RnaTargetClass

log = logging.getLogger(__name__)

SMRTNET_GITHUB = "https://github.com/Yuhan-Fei/SMRTnet.git"

_VALID_RNA = set("ACGU")
_DOTBRACKET_CHARS = set("().<>[]{}")

_TRAINING_FILE = "SMRTnet_data.txt"
_BENCHMARK_FILES = {
    "smrtnet_benchmark":      "SMRTnet_benchmark.txt",
    "smrtnet_naldb":          "SMRTnet_benchmark_NALDB.txt",
    "smrtnet_newpub":         "SMRTnet_benchmark_NewPub.txt",
    "smrtnet_rbind":          "SMRTnet_benchmark_RBIND.txt",
    "smrtnet_rsim":           "SMRTnet_benchmark_RSIM.txt",
    "smrtnet_smmrna":         "SMRTnet_benchmark_SMMRNA.txt",
}


def download(out_dir: Path) -> Path:
    """Sparse-clone the SMRTnet GitHub repo and copy data/ into ``out_dir``."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    repo_dir = out_dir.parent / "_smrtnet_clone"
    if not repo_dir.exists():
        log.info("Cloning %s → %s", SMRTNET_GITHUB, repo_dir)
        subprocess.run(
            ["git", "clone", "--depth", "1", SMRTNET_GITHUB, str(repo_dir)],
            check=True,
        )
    src = repo_dir / "data"
    for f in src.glob("*.txt"):
        (out_dir / f.name).write_bytes(f.read_bytes())
    return out_dir


# ---- helpers ----

def _clean_rna(s: str) -> str | None:
    if not isinstance(s, str):
        return None
    s = s.strip().upper().replace("T", "U")
    if not s or set(s) - _VALID_RNA:
        return None
    return s


def _classify_rna_target(seq: str, source_tag: str) -> RnaTargetClass:
    """SMRTnet doesn't always tell us the family, so we return UNKNOWN by
    default. The benchmark splits do come from labeled databases (R-BIND, etc.)
    where families are known, but those labels are encoded in the source tag.
    """

    return RnaTargetClass.OTHER


def _smi(smi: str) -> tuple[str, str] | None:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToSmiles(m), rd_inchi.MolToInchiKey(m)


def _is_dotbracket(s: str) -> bool:
    if not isinstance(s, str) or len(s) == 0:
        return False
    return set(s) <= _DOTBRACKET_CHARS


# ---- main loaders ----

def _load_one(path: Path, *, has_id_col: bool, source_tag: str) -> pd.DataFrame:
    """Generic 4 / 5-column SMRTnet TSV parser."""

    log.info("Reading %s (source=%s)", path.name, source_tag)
    rows: list[dict] = []
    n_skipped_smi = 0
    n_skipped_rna = 0
    n_skipped_db = 0
    with path.open("r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            if has_id_col:
                if len(parts) < 5:
                    continue
                _, smi, rna, dot, label = parts[:5]
            else:
                if len(parts) < 4:
                    continue
                smi, rna, dot, label = parts[:4]
            seq = _clean_rna(rna)
            if seq is None:
                n_skipped_rna += 1
                continue
            if not _is_dotbracket(dot) or len(dot) != len(seq):
                n_skipped_db += 1
                continue
            c = _smi(smi)
            if c is None:
                n_skipped_smi += 1
                continue
            canon, ikey = c
            try:
                lab = int(float(label))
            except ValueError:
                continue
            rows.append(
                dict(
                    record_id=f"{source_tag}_{i:06d}",
                    source=source_tag,
                    source_version="2026-01",
                    citation="https://doi.org/10.1038/s41587-025-02942-z",
                    rna_sequence=seq,
                    rna_class=_classify_rna_target(seq, source_tag).value,
                    rna_length=len(seq),
                    rna_secondary_structure=dot,
                    ligand_smiles=canon,
                    ligand_inchikey=ikey,
                    assay_kind=AssayKind.UNKNOWN.value,
                    activity_value=None,
                    activity_units=None,
                    is_active=bool(lab),
                )
            )
    log.info(
        "  %s: parsed %d rows | skipped: bad_smi=%d bad_rna=%d bad_db=%d",
        path.name, len(rows), n_skipped_smi, n_skipped_rna, n_skipped_db,
    )
    return pd.DataFrame(rows)


def load(raw_dir: Path) -> dict[str, pd.DataFrame]:
    """Load every SMRTnet TSV into a dict ``{source_tag: dataframe}``."""

    raw_dir = Path(raw_dir)
    out: dict[str, pd.DataFrame] = {}
    if (raw_dir / _TRAINING_FILE).exists():
        out["smrtnet_train"] = _load_one(
            raw_dir / _TRAINING_FILE,
            has_id_col=False,
            source_tag="smrtnet_train",
        )
    for tag, fname in _BENCHMARK_FILES.items():
        path = raw_dir / fname
        if not path.exists():
            log.warning("missing %s", path)
            continue
        out[tag] = _load_one(path, has_id_col=True, source_tag=tag)
    return out


def to_records(df: pd.DataFrame):
    for _, row in df.iterrows():
        try:
            yield InteractionRecord(**{k: v for k, v in row.items() if pd.notna(v)})
        except Exception as exc:  # noqa: BLE001
            log.warning("SMRTnet row failed validation: %s", exc)
