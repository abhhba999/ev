"""ROBIN loader.

ROBIN (Repository Of BInders to Nucleic acids) is a small molecule microarray
(SMM) screening library of 24,572 compounds tested against 36 nucleic acid
targets (25 RNA + 11 DNA). The Yazdani et al. 2022 paper provides public CSVs
for the per-target hits and the sequence table.

Code repo:  https://github.com/ky66/ROBIN
Paper:      Yazdani K, Jordan D, Yang M, Fullenkamp CR, Calabrese DR, Boer R,
            Hilimire T, Allen TEH, Khan RT, Schneekloth JS Jr.
            *Angew. Chem. Int. Ed.* (2023) 62, e202211358.
            doi: 10.1002/anie.202211358
Bulk data:  https://doi.org/10.6084/m9.figshare.20401974

We import three CSVs from the repo (small enough to commit-by-reference):

- ``SMM_Biomolecule_Hits.csv`` — per-compound binary "did this compound hit
  any RNA / any DNA / any G4" labels. Used to populate
  :class:`~boltzrna_diff.data.schema.BinaryBinderRecord`.
- ``SMM_Target_Hits.csv`` — per-(compound × target) binary hit matrix
  (24572 × 36). Used to produce per-target
  :class:`~boltzrna_diff.data.schema.InteractionRecord` rows.
- ``SMM_Sequence_Table_hit_rates.csv`` — per-target metadata: nucleic acid
  type, structural class, sequence, hit rate, selective hit rate.
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

from .schema import (
    AssayKind,
    BinaryBinderRecord,
    InteractionRecord,
    RnaTargetClass,
)

log = logging.getLogger(__name__)

ROBIN_GITHUB = "https://github.com/ky66/ROBIN.git"

# Map ROBIN's structural class strings to our coarse RnaTargetClass.
# ROBIN's "G-quadruplex" RNAs are tricky — we keep them as G_QUADRUPLEX.
_RNA_CLASS_MAP = {
    "G-quadruplex": RnaTargetClass.G_QUADRUPLEX,
    "Hairpin": RnaTargetClass.OTHER,  # default; we override below by name
    "Pseudoknot": RnaTargetClass.RIBOSWITCH,  # PreQ1, SAM_II, ZTP are riboswitches
    "Three-way junction": RnaTargetClass.RIBOSWITCH,  # TPP, Glutamine_RS
    "Triple helix": RnaTargetClass.LNCRNA,  # MALAT1, ENE_A9
}

# Per-target overrides because ROBIN's "Hairpin" bucket lumps pre-miRNAs and
# clinically interesting hairpins together. Be precise here.
_TARGET_NAME_OVERRIDES = {
    "Pre_miR_21": RnaTargetClass.PRE_MIRNA,
    "Pre_miR_17": RnaTargetClass.PRE_MIRNA,
    "Pre_miR_31": RnaTargetClass.PRE_MIRNA,
    "MALAT1": RnaTargetClass.LNCRNA,
    "TERRA": RnaTargetClass.G_QUADRUPLEX,
    "HIV_SL3": RnaTargetClass.VIRAL,
    "HBV": RnaTargetClass.VIRAL,
    "Zika3PrimeUTR": RnaTargetClass.VIRAL,
    "Zika_NS5": RnaTargetClass.VIRAL,
    "RRE2B": RnaTargetClass.VIRAL,
    "RRE2B_MeA": RnaTargetClass.VIRAL,
    "PreQ1": RnaTargetClass.RIBOSWITCH,
    "SAM_ll": RnaTargetClass.RIBOSWITCH,
    "TPP": RnaTargetClass.RIBOSWITCH,
    "ZTP": RnaTargetClass.RIBOSWITCH,
    "ENE_A9": RnaTargetClass.LNCRNA,
}


def download(out_dir: Path) -> Path:
    """Sparse-clone ROBIN's SMM_full_results into ``out_dir``."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    repo_dir = out_dir.parent / "_robin_clone"
    if not repo_dir.exists():
        log.info("Cloning %s → %s", ROBIN_GITHUB, repo_dir)
        subprocess.run(
            ["git", "clone", "--depth", "1", ROBIN_GITHUB, str(repo_dir)],
            check=True,
        )
    src = repo_dir / "SMM_full_results"
    for f in src.glob("*.csv"):
        (out_dir / f.name).write_bytes(f.read_bytes())
    return out_dir


# ---- helpers ----

_VALID_RNA = set("ACGU")


def _clean_rna_sequence(s: str) -> str | None:
    """ROBIN sometimes uses '(m6A)' modified bases or 'and' joined fragments.

    For the v0 schema we strip modifications back to the natural base and
    drop entries with split-fragment sequences.
    """

    if not isinstance(s, str):
        return None
    s = s.upper().replace("T", "U").strip()
    if " AND " in s:
        return None  # split-fragment, skip
    s = re.sub(r"\(M6A\)", "A", s)
    s = re.sub(r"\(M\d?[ACGU]\)", lambda m: m.group(0)[-2], s)  # generic (m1A) -> A
    if not set(s) <= _VALID_RNA:
        return None
    return s


def _smi_canon_inchikey(smi: str) -> tuple[str, str] | None:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToSmiles(m), rd_inchi.MolToInchiKey(m)


# ---- main loaders ----

def load_binary(raw_dir: Path) -> pd.DataFrame:
    """Build per-compound :class:`BinaryBinderRecord` rows."""

    raw_dir = Path(raw_dir)
    df = pd.read_csv(raw_dir / "SMM_Biomolecule_Hits.csv")
    rows: list[dict] = []
    for _, r in tqdm(df.iterrows(), total=len(df), desc="ROBIN binary"):
        smi = str(r["Smile"]).strip()
        c = _smi_canon_inchikey(smi)
        if c is None:
            continue
        smi, ikey = c
        rows.append(
            dict(
                record_id=f"robin_bin_{r['Name']}",
                source="robin",
                source_version="2022-08",
                citation="https://doi.org/10.1002/anie.202211358",
                ligand_smiles=smi,
                ligand_inchikey=ikey,
                is_rna_binder=bool(int(r["Hit_RNA"])),
                is_dna_binder=bool(int(r["Hit_DNA"])),
                assay_kind=AssayKind.MICROARRAY.value,
            )
        )
    return pd.DataFrame(rows)


def load_per_target(raw_dir: Path, rna_only: bool = True) -> pd.DataFrame:
    """Build per-target :class:`InteractionRecord` rows from ROBIN.

    If ``rna_only`` is True (default), DNA targets are dropped.
    """

    raw_dir = Path(raw_dir)
    hits = pd.read_csv(raw_dir / "SMM_Target_Hits.csv")
    seq = pd.read_csv(raw_dir / "SMM_Sequence_Table_hit_rates.csv")
    seq.columns = [c.strip() for c in seq.columns]

    target_meta: dict[str, dict] = {}
    for _, r in seq.iterrows():
        name = str(r["Nucleic Acid Target"]).strip()
        biomol = str(r["Biomolecule Type"]).strip().upper()
        if rna_only and biomol != "RNA":
            continue
        rna = _clean_rna_sequence(r["Sequence (5' to 3')"])
        if rna is None:
            continue
        struct = str(r["Structure Type"]).strip()
        cls = _TARGET_NAME_OVERRIDES.get(name, _RNA_CLASS_MAP.get(struct, RnaTargetClass.OTHER))
        target_meta[name] = dict(
            sequence=rna,
            length=len(rna),
            structure_type=struct,
            rna_class=cls.value,
        )
    log.info("ROBIN: %d RNA targets with usable sequences", len(target_meta))

    # The hit columns are named "<TargetName>_hit". Build a column→target map.
    target_columns = {}
    for col in hits.columns:
        if col.endswith("_hit"):
            target = col[:-len("_hit")]
            if target in target_meta:
                target_columns[col] = target

    rows: list[dict] = []
    for _, r in tqdm(hits.iterrows(), total=len(hits), desc="ROBIN per-target"):
        smi = str(r["Smile"]).strip()
        c = _smi_canon_inchikey(smi)
        if c is None:
            continue
        smi, ikey = c
        cmpd = str(r["Name"]).strip()
        for col, target_name in target_columns.items():
            val = r[col]
            if pd.isna(val):
                continue
            tm = target_meta[target_name]
            rows.append(
                dict(
                    record_id=f"robin_{cmpd}_{target_name}",
                    source="robin",
                    source_version="2022-08",
                    citation="https://doi.org/10.1002/anie.202211358",
                    rna_sequence=tm["sequence"],
                    rna_class=tm["rna_class"],
                    rna_length=tm["length"],
                    rna_secondary_structure=None,
                    ligand_smiles=smi,
                    ligand_inchikey=ikey,
                    assay_kind=AssayKind.MICROARRAY.value,
                    activity_value=None,
                    activity_units=None,
                    is_active=bool(int(val)),
                )
            )
    df = pd.DataFrame(rows)
    log.info(
        "ROBIN per-target: %d rows | %d unique compounds | %d targets | actives=%d",
        len(df),
        df["ligand_inchikey"].nunique() if len(df) else 0,
        df["rna_sequence"].nunique() if len(df) else 0,
        int(df["is_active"].sum()) if len(df) else 0,
    )
    return df


def to_binary_records(df: pd.DataFrame):
    for _, row in df.iterrows():
        try:
            yield BinaryBinderRecord(**{k: v for k, v in row.items() if pd.notna(v)})
        except Exception as exc:  # noqa: BLE001
            log.warning("ROBIN binary row failed validation: %s", exc)


def to_interaction_records(df: pd.DataFrame):
    for _, row in df.iterrows():
        payload = {k: v for k, v in row.items() if pd.notna(v)}
        try:
            yield InteractionRecord(**payload)
        except Exception as exc:  # noqa: BLE001
            log.warning("ROBIN interaction row failed validation: %s", exc)
