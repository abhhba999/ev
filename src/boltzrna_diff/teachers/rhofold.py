"""RhoFold+ wrapper for RNA-only ensemble sampling.

RhoFold+ (Shen et al., *Nature Methods* 2024, doi:10.1038/s41592-024-02487-0)
predicts RNA tertiary structures and supports stochastic ensemble sampling
via dropout-at-inference. It is our primary source of "conformation
ensembles" used to teach the student model that RNA targets are
conformationally heterogeneous.

We invoke the upstream ``rhofold`` CLI (https://github.com/ml4bio/RhoFold).
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path

from .types import EnsembleSample, TeacherInputs, TeacherOutputs

log = logging.getLogger(__name__)


def _detect_rhofold_binary() -> str | None:
    return shutil.which("rhofold")


def run(
    inputs: TeacherInputs,
    work_dir: Path,
    *,
    num_samples: int = 8,
    binary: str | None = None,
    extra_args: list[str] | None = None,
) -> TeacherOutputs:
    """Run RhoFold+ ensemble sampling.

    Returns one :class:`EnsembleSample` per dropout replica; the structures
    can then be fed to Boltz / AF3 for ligand-conditioned cofolding.
    """

    binary = binary or _detect_rhofold_binary()
    if binary is None:
        raise RuntimeError(
            "rhofold CLI not found. Install with `pip install rhofold` on a "
            "GPU node with CUDA 11.8+. See docs/cluster.md for SLURM template."
        )

    out_dir = Path(work_dir) / f"{inputs.target_id}_rhofold"
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / f"{inputs.target_id}.fasta"
    fasta.write_text(f">{inputs.target_id}\n{inputs.rna_sequence}\n")

    cmd = [
        binary,
        "predict",
        "--input_fasta",
        str(fasta),
        "--output_dir",
        str(out_dir),
        "--num_samples",
        str(num_samples),
        "--seed",
        str(inputs.seed),
    ]
    if extra_args:
        cmd += list(extra_args)
    log.info("$ %s", " ".join(cmd))
    subprocess.run(cmd, check=True, env={**os.environ})

    samples = []
    for i, p in enumerate(sorted(out_dir.glob("*.pdb"))):
        samples.append(
            EnsembleSample(
                target_id=inputs.target_id,
                sample_id=i,
                structure_path=p,
            )
        )
    return TeacherOutputs(
        target_id=inputs.target_id,
        teacher_name="rhofold+",
        samples=samples,
        output_dir=out_dir,
    )


SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=rhofold-{target_id}
#SBATCH --gres=gpu:a100:1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/rhofold_{target_id}_%j.out

set -euo pipefail
source ~/.bashrc
conda activate rhofold

uv run python -m boltzrna_diff.teachers.rhofold_run \\
    --target {target_id} \\
    --rna {rna_sequence} \\
    --num-samples {num_samples} \\
    --out {out_dir}
"""
