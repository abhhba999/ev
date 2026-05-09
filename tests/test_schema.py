"""Smoke tests for the unified RNA-rSM-Bench schema."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from boltzrna_diff.data.schema import (
    AssayKind,
    BinaryBinderRecord,
    InteractionRecord,
    RibotacRecord,
    RnaTargetClass,
    StructuralMethod,
    StructureRecord,
)


def test_structure_record_canonicalizes_smiles_and_uppercases_rna():
    r = StructureRecord(
        record_id="hariboss_1FMN_FMN",
        source="hariboss",
        source_version="2025-04",
        rna_sequence="GgCaUu".replace("T", "U"),
        rna_class=RnaTargetClass.RIBOSWITCH,
        rna_length=6,
        pdb_id="1FMN",
        structural_method=StructuralMethod.NMR,
        ligand_smiles="C1=CC=CC=C1",  # benzene, non-canonical
        ligand_inchikey="UHOVQNZJYSORNB-UHFFFAOYSA-N",
    )
    assert r.rna_sequence == "GGCAUU"
    assert r.ligand_smiles == "c1ccccc1"


def test_structure_record_rejects_bad_rna_alphabet():
    with pytest.raises(ValidationError):
        StructureRecord(
            record_id="bad",
            source="hariboss",
            source_version="2025-04",
            rna_sequence="GGGXYZ",
            rna_class=RnaTargetClass.OTHER,
            rna_length=6,
            ligand_smiles="CCO",
            ligand_inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        )


def test_structure_record_rejects_invalid_smiles():
    with pytest.raises(ValidationError):
        StructureRecord(
            record_id="bad",
            source="hariboss",
            source_version="2025-04",
            rna_sequence="GGCAUU",
            rna_class=RnaTargetClass.OTHER,
            rna_length=6,
            ligand_smiles="not a smiles@@@",
            ligand_inchikey="X",
        )


def test_interaction_record_round_trip():
    r = InteractionRecord(
        record_id="rbind_001",
        source="r_bind",
        source_version="2.0",
        rna_sequence="GGCAUU",
        rna_class=RnaTargetClass.PRE_MIRNA,
        rna_length=6,
        ligand_smiles="CCO",
        ligand_inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        assay_kind=AssayKind.ITC,
        activity_value=120.0,
        activity_units="nM",
        is_active=True,
    )
    j = r.model_dump_json()
    r2 = InteractionRecord.model_validate_json(j)
    assert r2.activity_value == 120.0


def test_binary_binder_record():
    r = BinaryBinderRecord(
        record_id="robin_001",
        source="robin",
        source_version="2024-12",
        ligand_smiles="CCO",
        ligand_inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        is_rna_binder=True,
        is_dna_binder=False,
    )
    assert r.assay_kind == AssayKind.MICROARRAY


def test_ribotac_record_canonicalizes_smiles():
    r = RibotacRecord(
        record_id="ribotac_costales_2020_1",
        source="ribotac_lit_v0",
        source_version="2026-04",
        rna_target_sequence="UAGCUUAUCAGACUGAUGUUGA",
        rna_class=RnaTargetClass.PRE_MIRNA,
        rna_target_name="pre-miR-21",
        smiles="C1=CC=CC=C1",
        inchikey="UHOVQNZJYSORNB-UHFFFAOYSA-N",
        dc50_nM=200.0,
        dmax_percent=80.0,
        cell_line="MDA-MB-231",
    )
    assert r.smiles == "c1ccccc1"
