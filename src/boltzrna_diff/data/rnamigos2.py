"""RNAmigos2 loader.

RNAmigos2 (Carvajal-Patino et al., *Nature Communications* 2025,
doi:10.1038/s41467-025-57852-0) is the first experimentally validated
structure-based deep-learning RNA virtual screening tool. The repository
ships a pre-packaged ``rnamigos2_data.tar.gz`` containing ~1,750 RNA pockets
(JSON graph form), 1.3 M binary native-vs-decoy labels, and 1.3 M docking
scores, all split into TRAIN / VALIDATION / TEST.

GitHub: https://github.com/cgoliver/rnamigos2

Files of interest (inside the tarball):

- ``csvs/binary_data.csv`` — 1,315,816 rows × ``[PDB_ID_POCKET, LIGAND_SMILES, IS_NATIVE, SPLIT]``
- ``csvs/docking_data.csv`` — same shape but with ``INTER`` (docking score) and ``normalized_values``
- ``csvs/fp_data.csv`` — 1,543 rows of fingerprint targets
- ``json_pockets_expanded/<pocket_id>.json`` — 1,749 NetworkX-style graphs
  with ``nodes = [{nt_code, id, in_pocket}, ...]`` so that we can recover
  the RNA sequence + the ≤6 Å pocket residue mask

Plus ``data/ROBIN.csv`` (in the bare repo root, not the tarball) — a small
9-row mapping ROBIN_ID → reference PDB / ligand for the ROBIN benchmark.

This loader produces ``InteractionRecord``-shaped rows (one per
(pocket × ligand × split-set) tuple). The ``rna_sequence`` is reconstructed
from the corresponding pocket JSON.
"""

from __future__ import annotations

import json
import logging
import re
import subprocess
import tarfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from tqdm import tqdm

from .schema import AssayKind, InteractionRecord, RnaTargetClass

log = logging.getLogger(__name__)

RNAMIGOS_GITHUB = "https://github.com/cgoliver/rnamigos2.git"

_VALID_RNA = set("ACGU")
_RES_ID_RE = re.compile(r"^([0-9a-zA-Z]+)\.([A-Za-z0-9]+)\.(-?\d+)$")


@dataclass
class PocketInfo:
    """Parsed RNA pocket from json_pockets_expanded/*.json."""

    pocket_id: str
    pdb_id: str
    chain: str
    sequence: str
    pocket_residues: list[int]  # 0-based positions inside the contiguous chain


def download(out_dir: Path) -> Path:
    """Sparse-clone the RNAmigos2 GitHub repo and copy data/ into ``out_dir``.

    The 57 MB tarball is included; we expand only ``csvs/`` and the
    ``json_pockets_expanded/`` directory.
    """

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    repo_dir = out_dir.parent / "_rnamigos2_clone"
    if not repo_dir.exists():
        log.info("Cloning %s → %s", RNAMIGOS_GITHUB, repo_dir)
        subprocess.run(
            ["git", "clone", "--depth", "1", RNAMIGOS_GITHUB, str(repo_dir)],
            check=True,
        )

    # Copy ROBIN.csv (small mapping table)
    src_robin = repo_dir / "data" / "ROBIN.csv"
    if src_robin.exists():
        (out_dir / "ROBIN.csv").write_bytes(src_robin.read_bytes())

    # Extract csvs/ and json_pockets_expanded/ from the tarball
    tarball = repo_dir / "data" / "rnamigos2_data.tar.gz"
    if tarball.exists() and not (out_dir / "csvs").exists():
        log.info("Extracting tarball → %s", out_dir)
        with tarfile.open(tarball, "r:gz") as tar:
            tar.extractall(out_dir)
    return out_dir


# ---- pocket parsing ----

def _parse_node_id(node_id: str) -> tuple[str, str, int] | None:
    m = _RES_ID_RE.match(node_id)
    if m is None:
        return None
    pdb_id, chain, pos = m.groups()
    return pdb_id, chain, int(pos)


def parse_pocket_json(path: Path) -> PocketInfo | None:
    with path.open("r", encoding="utf-8") as f:
        graph = json.load(f)
    nodes = graph.get("nodes", [])
    parsed: list[tuple[int, str, bool, str, str]] = []  # pos, chain, in_pocket, nt, pdb
    for n in nodes:
        nid = n.get("id", "")
        nt = n.get("nt_code", "").upper().replace("T", "U")
        ip = bool(n.get("in_pocket", False))
        info = _parse_node_id(nid)
        if info is None or nt not in _VALID_RNA:
            continue
        pdb_id, chain, pos = info
        parsed.append((pos, chain, ip, nt, pdb_id))
    if not parsed:
        return None
    # Pick the chain with most residues (pockets sometimes span chains)
    by_chain: dict[str, list[tuple[int, str, bool, str, str]]] = {}
    for row in parsed:
        by_chain.setdefault(row[1], []).append(row)
    main_chain = max(by_chain, key=lambda c: len(by_chain[c]))
    rows = sorted(by_chain[main_chain])
    sequence = "".join(r[3] for r in rows)
    pos_to_idx = {r[0]: i for i, r in enumerate(rows)}
    pocket_residues = [pos_to_idx[r[0]] for r in rows if r[2]]
    pdb_id = rows[0][4]
    return PocketInfo(
        pocket_id=path.stem,
        pdb_id=pdb_id,
        chain=main_chain,
        sequence=sequence,
        pocket_residues=pocket_residues,
    )


def load_pockets(raw_dir: Path) -> dict[str, PocketInfo]:
    raw_dir = Path(raw_dir)
    pdir = raw_dir / "json_pockets_expanded"
    if not pdir.exists():
        raise FileNotFoundError(f"missing {pdir} — run download() first")
    pockets: dict[str, PocketInfo] = {}
    for path in tqdm(sorted(pdir.glob("*.json")), desc="parsing pockets"):
        info = parse_pocket_json(path)
        if info is not None:
            pockets[info.pocket_id] = info
    log.info("Parsed %d pockets", len(pockets))
    return pockets


# ---- main loaders ----

def _smi(smi: str) -> tuple[str, str] | None:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToSmiles(m), rd_inchi.MolToInchiKey(m)


def _build_rows(
    df: pd.DataFrame,
    pockets: dict[str, PocketInfo],
    *,
    source_tag: str,
    activity_col: str | None,
    is_active_col: str | None,
    activity_units: str | None = None,
) -> pd.DataFrame:
    out: list[dict] = []
    n_no_pocket = 0
    n_bad_smi = 0
    for i, row in tqdm(df.iterrows(), total=len(df), desc=f"build {source_tag}"):
        pid = row["PDB_ID_POCKET"]
        if pid not in pockets:
            n_no_pocket += 1
            continue
        pkt = pockets[pid]
        c = _smi(row["LIGAND_SMILES"])
        if c is None:
            n_bad_smi += 1
            continue
        canon, ikey = c
        rec = {
            "record_id": f"{source_tag}_{i}",
            "source": source_tag,
            "source_version": "2025-03",
            "citation": "https://doi.org/10.1038/s41467-025-57852-0",
            "rna_sequence": pkt.sequence,
            "rna_class": RnaTargetClass.OTHER.value,
            "rna_length": len(pkt.sequence),
            "rna_secondary_structure": None,
            "ligand_smiles": canon,
            "ligand_inchikey": ikey,
            "assay_kind": AssayKind.UNKNOWN.value,
            "activity_value": float(row[activity_col]) if activity_col else None,
            "activity_units": activity_units,
            "is_active": bool(row[is_active_col]) if is_active_col else None,
            "pdb_link": pkt.pdb_id,
        }
        # Stash the official RNAmigos2 SPLIT column as a side channel
        if "SPLIT" in row.index:
            rec["_rnamigos_split"] = row["SPLIT"]
        out.append(rec)
    log.info(
        "  %s: built %d rows | skipped: no_pocket=%d bad_smi=%d",
        source_tag, len(out), n_no_pocket, n_bad_smi,
    )
    return pd.DataFrame(out)


def load_binary(raw_dir: Path, pockets: dict[str, PocketInfo] | None = None) -> pd.DataFrame:
    """Native-vs-decoy classification (1.3 M rows)."""

    raw_dir = Path(raw_dir)
    if pockets is None:
        pockets = load_pockets(raw_dir)
    df = pd.read_csv(raw_dir / "csvs" / "binary_data.csv", index_col=0)
    return _build_rows(
        df,
        pockets,
        source_tag="rnamigos2_binary",
        activity_col=None,
        is_active_col="IS_NATIVE",
    )


def load_docking(raw_dir: Path, pockets: dict[str, PocketInfo] | None = None) -> pd.DataFrame:
    """Docking-score regression labels (1.3 M rows, INTER in kcal/mol)."""

    raw_dir = Path(raw_dir)
    if pockets is None:
        pockets = load_pockets(raw_dir)
    df = pd.read_csv(raw_dir / "csvs" / "docking_data.csv", index_col=0)
    return _build_rows(
        df,
        pockets,
        source_tag="rnamigos2_docking",
        activity_col="INTER",
        is_active_col=None,
        activity_units="kcal/mol",
    )


def load_robin_targets(raw_dir: Path) -> pd.DataFrame:
    """Read the small ROBIN.csv mapping (9 rows: ROBIN_ID → PDB / ligand)."""

    return pd.read_csv(Path(raw_dir) / "ROBIN.csv")


def to_records(df: pd.DataFrame):
    for _, row in df.iterrows():
        try:
            payload = {
                k: v for k, v in row.items()
                if pd.notna(v) and not k.startswith("_")
            }
            yield InteractionRecord(**payload)
        except Exception as exc:  # noqa: BLE001
            log.warning("RNAmigos2 row failed validation: %s", exc)
