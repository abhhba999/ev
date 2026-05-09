"""Smoke tests for teacher-model wrappers (CPU-only path)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from boltzrna_diff.teachers.boltz import build_yaml_input, render_slurm_script
from boltzrna_diff.teachers.types import TeacherInputs


def test_boltz_yaml_with_smiles(tmp_path: Path):
    inputs = TeacherInputs(
        target_id="pre_miR_21",
        rna_sequence="GGGUUGACUGUUGAAUCUCAUGGCAACCC",
        ligand_smiles="O=C(O)Cc1ccc(N)cc1",
    )
    p = build_yaml_input(inputs, tmp_path)
    payload = json.loads(p.read_text())
    seqs = payload["sequences"]
    assert seqs[0]["rna"]["sequence"] == "GGGUUGACUGUUGAAUCUCAUGGCAACCC"
    assert seqs[1]["ligand"]["smiles"] == "O=C(O)Cc1ccc(N)cc1"


def test_boltz_yaml_with_sdf(tmp_path: Path):
    sdf = tmp_path / "lig.sdf"
    sdf.write_text("dummy")
    inputs = TeacherInputs(
        target_id="fmn",
        rna_sequence="GGAUAUGAGGCG",
        ligand_sdf_path=sdf,
    )
    p = build_yaml_input(inputs, tmp_path)
    payload = json.loads(p.read_text())
    assert payload["sequences"][1]["ligand"]["sdf"].endswith("lig.sdf")


def test_boltz_run_raises_without_binary(monkeypatch, tmp_path: Path):
    from boltzrna_diff.teachers import boltz

    monkeypatch.setattr(boltz, "_detect_boltz_binary", lambda: None)
    inputs = TeacherInputs(target_id="x", rna_sequence="GCAUG")
    with pytest.raises(RuntimeError, match="boltz CLI not found"):
        boltz.run(inputs, tmp_path)


def test_boltz_slurm_script_contains_target():
    s = render_slurm_script(
        target_id="hiv1_tar",
        rna_sequence="GGCAGAUCUGAGCCUGGGAGCUCUCUGCC",
        ligand_smiles="C1=CC=CC=C1",
        model="boltz1",
        num_samples=8,
    )
    assert "hiv1_tar" in s
    assert "boltz1" in s
    assert "8" in s


def test_rhofold_run_raises_without_binary(monkeypatch, tmp_path: Path):
    from boltzrna_diff.teachers import rhofold

    monkeypatch.setattr(rhofold, "_detect_rhofold_binary", lambda: None)
    inputs = TeacherInputs(target_id="x", rna_sequence="GCAUG")
    with pytest.raises(RuntimeError, match="rhofold CLI not found"):
        rhofold.run(inputs, tmp_path)
