"""Unified schema for the RNA-rSM-Bench v0.

Every downstream loader (HARIBOSS, ROBIN, R-BIND, Inforna, RNAmigos2-supp,
SMRTnet-supp, RIBOTAC literature) is responsible for emitting rows that
conform to one of these Pydantic models. We use Pydantic for validation
and pyarrow/polars for storage on disk (parquet).

There are three logical record types:

1. ``StructureRecord`` — an RNA + ligand 3D complex from PDB-derived sources
   (HARIBOSS) or computational complex datasets (RiboBind).
2. ``InteractionRecord`` — an RNA + ligand affinity / activity assay readout
   (R-BIND, RNAmigos2-supp, SMRTnet-supp, Inforna).
3. ``BinaryBinderRecord`` — an RNA-binder vs decoy classification entry
   (ROBIN-style).

Plus a separate ``RibotacRecord`` for bifunctional degraders that needs
extra metadata (RNase L recruiter, linker, DC50).
"""

from __future__ import annotations

from datetime import date
from enum import Enum
from typing import Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator
from rdkit import Chem


class RnaTargetClass(str, Enum):
    """Coarse RNA target taxonomy used for `family` OOD splits."""

    PRE_MIRNA = "pre_mirna"
    MIRNA = "mirna"
    LNCRNA = "lncrna"
    RIBOSWITCH = "riboswitch"
    APTAMER = "aptamer"
    IRES = "ires"
    PRE_MRNA_SS = "pre_mrna_splice_site"
    UTR_3 = "3_utr"
    UTR_5 = "5_utr"
    G_QUADRUPLEX = "g_quadruplex"
    REPEAT_EXPANSION = "repeat_expansion"
    VIRAL = "viral_rna"
    RIBOSOMAL = "ribosomal"
    OTHER = "other"


class AssayKind(str, Enum):
    """How was the activity / affinity measured."""

    ITC = "itc"
    NMR = "nmr"
    FRET = "fret"
    FLUORESCENCE_QUENCH = "fluorescence_quench"
    SPR = "spr"
    BLI = "bli"
    MS = "mass_spec"
    MICROARRAY = "microarray"
    DSF = "dsf"
    THERMAL_MELT = "thermal_melt"
    CELL_REPORTER = "cell_reporter"
    QPCR = "qpcr"
    SPLICING_ASSAY = "splicing_assay"
    DEGRADATION_DC50 = "degradation_dc50"
    UNKNOWN = "unknown"


class StructuralMethod(str, Enum):
    XRAY = "xray"
    NMR = "nmr"
    CRYOEM = "cryoem"
    AF3 = "af3_predicted"
    BOLTZ = "boltz_predicted"
    RHOFOLD = "rhofold_predicted"
    OTHER = "other"


class _CommonMixin(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    record_id: str = Field(..., description="Unique row id, e.g. 'hariboss_3F4G_LIG1'.")
    source: str = Field(..., description="Loader that emitted this row, e.g. 'hariboss'.")
    source_version: str = Field(..., description="Snapshot version, e.g. '2024-08-01'.")
    citation: Optional[str] = Field(None, description="DOI / URL / PubMed id of the source.")


class StructureRecord(_CommonMixin):
    """A 3D RNA-ligand complex."""

    rna_sequence: str = Field(..., description="Single-stranded RNA sequence (ACGU only).")
    rna_class: RnaTargetClass
    rna_length: int

    pdb_id: Optional[str] = Field(None, description="PDB accession, e.g. '6HAG'.")
    structural_method: StructuralMethod = StructuralMethod.OTHER
    resolution_angstrom: Optional[float] = None
    deposit_date: Optional[date] = None

    ligand_smiles: str = Field(..., description="Canonical SMILES of the bound ligand.")
    ligand_inchikey: str
    ligand_mw: Optional[float] = None

    pocket_residues: Optional[list[int]] = Field(
        default=None,
        description="0-based RNA residue indices forming the binding pocket (≤6 Å).",
    )
    structure_path: Optional[str] = Field(
        None, description="Relative path inside data/processed/ to the cleaned pdb/cif file."
    )

    @field_validator("rna_sequence")
    @classmethod
    def _check_rna(cls, v: str) -> str:
        v = v.upper().replace("T", "U")
        if not set(v) <= set("ACGU"):
            raise ValueError(f"non-ACGU character in rna_sequence: {v!r}")
        return v

    @field_validator("ligand_smiles")
    @classmethod
    def _canonicalize_smiles(cls, v: str) -> str:
        m = Chem.MolFromSmiles(v)
        if m is None:
            raise ValueError(f"invalid SMILES: {v!r}")
        return Chem.MolToSmiles(m)


class InteractionRecord(_CommonMixin):
    """A measured RNA + ligand affinity / activity readout."""

    rna_sequence: str
    rna_class: RnaTargetClass
    rna_length: int
    rna_secondary_structure: Optional[str] = Field(
        None, description="Dot-bracket notation if known."
    )

    ligand_smiles: str
    ligand_inchikey: str

    assay_kind: AssayKind
    activity_value: Optional[float] = Field(
        None,
        description="Numerical readout (Kd, EC50, IC50, DC50, fold-change, ...) in `activity_units`.",
    )
    activity_units: Optional[str] = Field(None, description="e.g. 'nM', 'uM', 'fold'.")
    is_active: Optional[bool] = Field(
        None, description="Active/inactive when the source reports a binary call."
    )

    pdb_link: Optional[str] = Field(
        None,
        description="If a parallel structural record exists in HARIBOSS, link it here.",
    )

    @field_validator("rna_sequence")
    @classmethod
    def _check_rna(cls, v: str) -> str:
        v = v.upper().replace("T", "U")
        if not set(v) <= set("ACGU.()"):
            raise ValueError(f"non-RNA character in rna_sequence: {v!r}")
        return v


class BinaryBinderRecord(_CommonMixin):
    """ROBIN-style binary binder/decoy table.

    Notable difference vs InteractionRecord: there is no specific RNA target —
    the data was generated by pooling many RNAs and asking 'does this molecule
    bind RNA at all?'.
    """

    ligand_smiles: str
    ligand_inchikey: str
    is_rna_binder: bool
    is_dna_binder: Optional[bool] = None
    assay_kind: AssayKind = AssayKind.MICROARRAY


class RibotacRecord(_CommonMixin):
    """Bifunctional RNA degrader (RIBOTAC).

    A RIBOTAC molecule = `rna_binder_smiles` + `linker_smiles` +
    `rnase_l_recruiter_smiles`. We always also store the *full* molecule
    (`smiles`) for round-trip checks.
    """

    rna_target_sequence: str
    rna_class: RnaTargetClass
    rna_target_name: str = Field(..., description="e.g. 'pre-miR-21', 'pre-miR-96'.")

    smiles: str = Field(..., description="Full RIBOTAC SMILES.")
    inchikey: str
    rna_binder_smiles: Optional[str] = Field(
        None,
        description="Sub-SMILES of the RNA-binding warhead, if dissectable.",
    )
    linker_smiles: Optional[str] = None
    rnase_l_recruiter_smiles: Optional[str] = None
    recruiter_kind: Optional[str] = Field(
        None,
        description="e.g. 'Kethoxal-recruiter', '2'-5'A4-mimic'.",
    )

    dc50_nM: Optional[float] = Field(None, description="Half-maximal degradation conc. (nM).")
    dmax_percent: Optional[float] = Field(None, description="Maximum % degradation observed.")
    cell_line: Optional[str] = None
    assay_kind: AssayKind = AssayKind.DEGRADATION_DC50

    @field_validator("smiles")
    @classmethod
    def _canon(cls, v: str) -> str:
        m = Chem.MolFromSmiles(v)
        if m is None:
            raise ValueError(f"invalid SMILES: {v!r}")
        return Chem.MolToSmiles(m)


__all__ = [
    "AssayKind",
    "BinaryBinderRecord",
    "InteractionRecord",
    "RibotacRecord",
    "RnaTargetClass",
    "StructuralMethod",
    "StructureRecord",
]
