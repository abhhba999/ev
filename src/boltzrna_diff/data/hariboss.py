"""HARIBOSS loader.

HARIBOSS (https://hariboss.pasteur.cloud/) is a curated database of RNA-small
molecule structures retrieved from the PDB. The website exposes two CSV
endpoints:

- ``/compounds/?format=csv``  → 326 unique ligands with physico-chemical props.
- ``/complexes/?format=csv``  → 1000 RNA-ligand complexes with metadata
  (PDB id, deposition date, experimental method, RNA chain sequences,
  ligand ids, pocket sizes, ...).

This loader joins the two CSVs into a list of ``StructureRecord``s.

Reference: Panei et al., Bioinformatics 38 (2022) 4185-4193.
DOI: https://doi.org/10.1093/bioinformatics/btac483
"""

from __future__ import annotations

import ast
import logging
from pathlib import Path
from typing import Iterator

import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from tqdm import tqdm

from .schema import RnaTargetClass, StructuralMethod, StructureRecord

log = logging.getLogger(__name__)

HARIBOSS_COMPOUNDS_URL = "https://hariboss.pasteur.cloud/compounds/?format=csv"
HARIBOSS_COMPLEXES_URL = "https://hariboss.pasteur.cloud/complexes/?format=csv"


# ---------- Download ----------

def download(out_dir: Path) -> tuple[Path, Path]:
    """Fetch the two HARIBOSS CSVs to ``out_dir`` and return their paths."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    compounds_path = out_dir / "compounds.csv"
    complexes_path = out_dir / "complexes.csv"
    headers = {"User-Agent": "Mozilla/5.0 (boltzrna-diff data loader)"}
    for url, path in [
        (HARIBOSS_COMPOUNDS_URL, compounds_path),
        (HARIBOSS_COMPLEXES_URL, complexes_path),
    ]:
        log.info("Downloading %s → %s", url, path)
        r = requests.get(url, headers=headers, timeout=120)
        r.raise_for_status()
        path.write_bytes(r.content)
    return compounds_path, complexes_path


# ---------- Helpers ----------

def _safe_literal(value: object) -> object:
    """HARIBOSS CSV cells embed Python-repr-encoded dicts/lists. Parse safely."""

    if not isinstance(value, str) or value.strip() == "":
        return None
    try:
        return ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return value


def _classify_rna(title: str, genus_field: object, sequence: str) -> RnaTargetClass:
    """Best-effort taxonomy from HARIBOSS free-text fields.

    HARIBOSS does not store an explicit RNA class column, so we use string
    heuristics on the title + (parsed) genus dict. This is intentionally coarse
    — the `family` OOD split will use this column.
    """

    t = (title or "").upper()
    if any(k in t for k in ("RIBOSWITCH", "FMN", "TPP", "GUANINE", "ADENINE", "SAM-")):
        return RnaTargetClass.RIBOSWITCH
    if "TAR" in t and "HIV" in t:
        return RnaTargetClass.VIRAL
    if "APTAMER" in t:
        return RnaTargetClass.APTAMER
    if any(k in t for k in ("MIRNA", "MIR-", "PRE-MIR", "MICRORNA")):
        return RnaTargetClass.PRE_MIRNA
    if any(k in t for k in ("RIBOSOM", "16S", "23S", "30S", "50S", "70S", "A-SITE")):
        return RnaTargetClass.RIBOSOMAL
    if any(k in t for k in ("IRES",)):
        return RnaTargetClass.IRES
    if any(k in t for k in ("REPEAT", "CUG", "CGG", "CAG", "G4C2")):
        return RnaTargetClass.REPEAT_EXPANSION
    if any(k in t for k in ("VIRAL", "HIV", "HEPATITIS", "CORONAVIRUS", "SARS", "DENGUE", "ZIKA")):
        return RnaTargetClass.VIRAL
    if any(k in t for k in ("G-QUADRUPLEX", "G-QUARTET", "G4")):
        return RnaTargetClass.G_QUADRUPLEX
    return RnaTargetClass.OTHER


def _classify_method(method: str | float | None) -> StructuralMethod:
    """HARIBOSS uses short codes: 'XR', 'EM', 'NMR'. Be permissive."""

    if not isinstance(method, str):
        return StructuralMethod.OTHER
    m = method.upper()
    if "XR" in m or "X-RAY" in m or m == "XRAY":
        return StructuralMethod.XRAY
    if "NMR" in m:
        return StructuralMethod.NMR
    if "EM" in m or "CRYO" in m:
        return StructuralMethod.CRYOEM
    return StructuralMethod.OTHER


def _ligand_inchikey(smiles: str) -> str | None:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return rd_inchi.MolToInchiKey(m)


# ---------- Main loader ----------

def load(raw_dir: Path) -> pd.DataFrame:
    """Load HARIBOSS into a flat dataframe of (complex × ligand) rows.

    Returns columns matching :class:`StructureRecord`.
    """

    raw_dir = Path(raw_dir)
    compounds = pd.read_csv(raw_dir / "compounds.csv")
    complexes = pd.read_csv(raw_dir / "complexes.csv")

    # Build SMILES lookup: ligand id (PDB chemical id) → row.
    cmp_lookup = compounds.set_index("id").to_dict(orient="index")

    rows: list[dict] = []
    for _, complex_row in tqdm(complexes.iterrows(), total=len(complexes), desc="HARIBOSS"):
        title = str(complex_row.get("title") or "")
        chains_raw = _safe_literal(complex_row.get("rna_chain_sequences"))
        if not isinstance(chains_raw, dict) or not chains_raw:
            continue
        # Concatenate chain sequences (HARIBOSS stores per-chain).
        seqs = [str(v).upper().replace("T", "U") for v in chains_raw.values() if v]
        seq = "".join(s for s in seqs if set(s) <= set("ACGU"))
        if not seq:
            continue
        rna_class = _classify_rna(title, complex_row.get("rna_genus"), seq)
        method = _classify_method(complex_row.get("experimental_method"))
        try:
            res = float(complex_row.get("experimental_resolution"))
        except (TypeError, ValueError):
            res = None
        deposit_date = complex_row.get("deposition_date")

        ligand_field = _safe_literal(complex_row.get("sm_ligand_ids"))
        ligand_pocket = _safe_literal(complex_row.get("sm_ligand_pocket_size")) or {}
        if not isinstance(ligand_field, list):
            continue

        for ligand_token in ligand_field:
            # ligand_token is e.g. "ARG_.:B/47:A". The molecular id (used in the
            # compounds table) is the prefix before the first underscore.
            ligand_id = ligand_token.split("_", 1)[0].strip()
            cmp_row = cmp_lookup.get(ligand_id)
            if not cmp_row:
                continue
            smiles = cmp_row.get("canonical_smiles")
            if not isinstance(smiles, str) or not smiles:
                continue
            inchikey = cmp_row.get("inchikey") or _ligand_inchikey(smiles)
            if inchikey is None:
                continue
            try:
                pocket_size = float(ligand_pocket.get(ligand_token)) if isinstance(ligand_pocket, dict) else None
            except (TypeError, ValueError):
                pocket_size = None
            rows.append(
                dict(
                    record_id=f"hariboss_{complex_row['id']}_{ligand_id}",
                    source="hariboss",
                    source_version="2025-04",
                    citation="https://doi.org/10.1093/bioinformatics/btac483",
                    rna_sequence=seq,
                    rna_class=rna_class.value,
                    rna_length=len(seq),
                    pdb_id=str(complex_row["id"]).upper(),
                    structural_method=method.value,
                    resolution_angstrom=res,
                    deposit_date=deposit_date if isinstance(deposit_date, str) else None,
                    ligand_smiles=smiles,
                    ligand_inchikey=str(inchikey),
                    ligand_mw=cmp_row.get("molecular_weight"),
                    pocket_size_a3=pocket_size,
                    structure_path=None,
                )
            )

    df = pd.DataFrame(rows)
    log.info("HARIBOSS produced %d (complex × ligand) rows from %d complexes / %d compounds",
             len(df), len(complexes), len(compounds))
    return df


def to_records(df: pd.DataFrame) -> Iterator[StructureRecord]:
    """Validate every row through the Pydantic schema."""

    drop_cols = {"pocket_size_a3"}  # extra column not in schema
    for _, row in df.iterrows():
        payload = {k: v for k, v in row.items() if k not in drop_cols and pd.notna(v)}
        try:
            yield StructureRecord(**payload)
        except Exception as exc:  # noqa: BLE001 — keep going on schema violations
            log.warning("HARIBOSS row %s failed validation: %s", row.get("record_id"), exc)
