"""Shared dataclasses for teacher-model wrappers."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class TeacherInputs:
    """Minimal inputs for a single (RNA target × ligand) inference call."""

    target_id: str
    rna_sequence: str
    ligand_smiles: Optional[str] = None
    ligand_sdf_path: Optional[Path] = None
    rna_msa_path: Optional[Path] = None
    pocket_residues: Optional[list[int]] = None
    seed: int = 0


@dataclass
class EnsembleSample:
    """One conformation produced by the teacher."""

    target_id: str
    sample_id: int
    structure_path: Path
    pae_path: Optional[Path] = None
    plddt_path: Optional[Path] = None
    interaction_score: Optional[float] = None  # for Boltz-2 affinity head
    extras: dict = field(default_factory=dict)


@dataclass
class TeacherOutputs:
    """All ensembles produced for a target by a single teacher."""

    target_id: str
    teacher_name: str
    samples: list[EnsembleSample]
    output_dir: Path
