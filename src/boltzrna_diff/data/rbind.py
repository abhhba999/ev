"""R-BIND 2.0 loader (manual-download).

R-BIND 2.0 (Donlic et al., *ACS Chem Biol* 2022, doi:10.1021/acschembio.2c00224)
is the curated database of bioactive RNA-targeting small molecules with
demonstrated cellular / animal activity. As of 2025, the website
(rbind.chem.duke.edu) is in transition (Hargrove Lab moved to Toronto) and the
ACS supplementary data is paywalled, so we cannot auto-download.

Workflow:
1. Manually download the ACS Chem Biol supplementary file
   ``cb2c00224_si_001.pdf`` and ``cb2c00224_si_002.xlsx`` from
   https://pubs.acs.org/doi/10.1021/acschembio.2c00224 (institutional access).
2. Place them in ``data/raw/rbind/``.
3. Run ``boltzrna-diff data build-rbind``.

The xlsx has columns roughly:
- Compound number / name
- SMILES (Canonical, ChemDraw)
- RNA target (free text)
- RNA secondary structure (free text)
- Activity (Kd / EC50 / IC50, often with units)
- Reference DOI

We map this to :class:`~boltzrna_diff.data.schema.InteractionRecord`. Since
the SI uses heterogeneous units, this loader normalizes everything to nM.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi

from .schema import AssayKind, InteractionRecord, RnaTargetClass

log = logging.getLogger(__name__)


_UNIT_TO_nM = {
    "pM": 1e-3,
    "nM": 1.0,
    "μM": 1e3,
    "uM": 1e3,
    "mM": 1e6,
    "M": 1e9,
}


def _parse_activity(text: str) -> tuple[float | None, str]:
    """Pull a numeric value + unit from free-text like '120 nM' or 'Kd = 1.2 µM'."""

    if not isinstance(text, str):
        return None, "unknown"
    m = re.search(r"([\d\.]+)\s*([numµμ]?[pnumμ]?M)", text)
    if m is None:
        return None, "unknown"
    value = float(m.group(1))
    unit = m.group(2).replace("μ", "u").replace("µ", "u")
    factor = _UNIT_TO_nM.get(unit, None)
    if factor is None:
        return None, unit
    return value * factor, "nM"


def _classify_rna_target(text: str) -> RnaTargetClass:
    if not isinstance(text, str):
        return RnaTargetClass.OTHER
    t = text.upper()
    if any(k in t for k in ("PRE-MIR", "PRE_MIR", "MIR-")):
        return RnaTargetClass.PRE_MIRNA
    if "RIBOSWITCH" in t or any(k in t for k in ("FMN", "TPP", "SAM", "PREQ1", "GLN", "ZTP")):
        return RnaTargetClass.RIBOSWITCH
    if "TAR" in t or "HIV" in t or "HCV" in t or "FRAMESHIFT" in t:
        return RnaTargetClass.VIRAL
    if "MALAT" in t or "NEAT" in t or "LINCRNA" in t or "LNCRNA" in t:
        return RnaTargetClass.LNCRNA
    if "G-QUAD" in t or "G4" in t:
        return RnaTargetClass.G_QUADRUPLEX
    if "REPEAT" in t or any(k in t for k in ("CUG", "CGG", "CAG", "G4C2")):
        return RnaTargetClass.REPEAT_EXPANSION
    if "5'SS" in t or "5'-SS" in t or "SPLIC" in t:
        return RnaTargetClass.PRE_MRNA_SS
    if "RIBOSOMAL" in t or "16S" in t or "23S" in t or "A-SITE" in t:
        return RnaTargetClass.RIBOSOMAL
    return RnaTargetClass.OTHER


def load(raw_dir: Path) -> pd.DataFrame:
    """Parse R-BIND 2.0 SI xlsx into our schema.

    Expects ``raw_dir/cb2c00224_si_002.xlsx`` (or any xlsx). Tolerates column-
    name variants by lowercasing and matching keywords.
    """

    raw_dir = Path(raw_dir)
    candidates = list(raw_dir.glob("*.xlsx"))
    if not candidates:
        raise FileNotFoundError(
            f"No R-BIND xlsx found in {raw_dir}. See module docstring for manual download."
        )
    xlsx = candidates[0]
    log.info("Reading R-BIND from %s", xlsx)
    df = pd.read_excel(xlsx)

    cols = {c.lower(): c for c in df.columns}
    smi_col = next((cols[k] for k in cols if "smile" in k), None)
    seq_col = next((cols[k] for k in cols if "sequence" in k or "rna" in k and "target" not in k), None)
    target_col = next((cols[k] for k in cols if "target" in k or "biomol" in k), None)
    activity_col = next((cols[k] for k in cols if "activity" in k or "kd" in k or "ic50" in k or "ec50" in k), None)
    ref_col = next((cols[k] for k in cols if "doi" in k or "ref" in k), None)

    if smi_col is None:
        raise ValueError(f"Could not find SMILES column in {xlsx} (have {list(df.columns)})")

    rows: list[dict] = []
    for i, r in df.iterrows():
        smi = str(r[smi_col]).strip()
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        canon = Chem.MolToSmiles(m)
        ikey = rd_inchi.MolToInchiKey(m)
        target_text = str(r.get(target_col, "") or "")
        seq = str(r.get(seq_col, "") or "")
        seq_clean = "".join(ch for ch in seq.upper().replace("T", "U") if ch in "ACGU")
        if not seq_clean:
            continue
        cls = _classify_rna_target(target_text + " " + seq)
        act_value, act_unit = _parse_activity(str(r.get(activity_col, "") or ""))
        rows.append(
            dict(
                record_id=f"rbind_{i+1:04d}",
                source="r_bind",
                source_version="2.0",
                citation=str(r.get(ref_col, "") or "https://doi.org/10.1021/acschembio.2c00224"),
                rna_sequence=seq_clean,
                rna_class=cls.value,
                rna_length=len(seq_clean),
                rna_secondary_structure=None,
                ligand_smiles=canon,
                ligand_inchikey=ikey,
                assay_kind=AssayKind.UNKNOWN.value,
                activity_value=act_value,
                activity_units=act_unit,
                is_active=True,
            )
        )
    log.info("R-BIND loaded %d rows", len(rows))
    return pd.DataFrame(rows)


def to_records(df: pd.DataFrame):
    for _, row in df.iterrows():
        try:
            yield InteractionRecord(**{k: v for k, v in row.items() if pd.notna(v)})
        except Exception as exc:  # noqa: BLE001
            log.warning("R-BIND row failed validation: %s", exc)
