"""Boltz-1 / Boltz-2 cofolding wrapper.

This module is imported on the CPU workstation for *config + IO* only; the
actual neural-network inference must run on the GPU cluster (see
``docs/cluster.md``).

Usage on a CPU box (e.g. for unit tests / DAG testing)::

    from boltzrna_diff.teachers.boltz import build_yaml_input
    yaml_path = build_yaml_input(inputs, work_dir)

Usage on a GPU node (after ``conda activate boltz``)::

    from boltzrna_diff.teachers.boltz import run
    out = run(inputs, work_dir, model="boltz1", num_samples=16)

We deliberately wrap Boltz via the ``boltz predict`` CLI rather than the
Python API, because the CLI handles MSA fetching, model checkpoints, and
diffusion sampling end-to-end. This is the same pattern used by Carvajal-
Patiño et al. for RNAmigos2 inference.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Literal

from .types import EnsembleSample, TeacherInputs, TeacherOutputs

log = logging.getLogger(__name__)

ModelName = Literal["boltz1", "boltz2"]

# ---- input building ----

def build_yaml_input(inputs: TeacherInputs, work_dir: Path) -> Path:
    """Write the Boltz YAML input file for one target.

    Schema reference:
    https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
    """

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = work_dir / f"{inputs.target_id}.yaml"
    rna_block = {
        "rna": {
            "id": "A",
            "sequence": inputs.rna_sequence,
            "msa": str(inputs.rna_msa_path) if inputs.rna_msa_path else "empty",
        }
    }
    sequences = [rna_block]
    if inputs.ligand_smiles is not None:
        sequences.append({"ligand": {"id": "L", "smiles": inputs.ligand_smiles}})
    elif inputs.ligand_sdf_path is not None:
        sequences.append({"ligand": {"id": "L", "sdf": str(inputs.ligand_sdf_path)}})
    payload = {
        "version": 1,
        "sequences": sequences,
        "constraints": {},
    }
    # Boltz uses YAML; we emit JSON as a strict subset for stability across versions.
    yaml_path.write_text(json.dumps(payload, indent=2))
    return yaml_path


# ---- runner ----

def _detect_boltz_binary() -> str | None:
    return shutil.which("boltz")


def run(
    inputs: TeacherInputs,
    work_dir: Path,
    *,
    model: ModelName = "boltz1",
    num_samples: int = 16,
    use_msa_server: bool = True,
    binary: str | None = None,
    extra_args: list[str] | None = None,
) -> TeacherOutputs:
    """Run Boltz-1/2 prediction on a GPU node.

    Returns
    -------
    TeacherOutputs containing one :class:`EnsembleSample` per requested
    diffusion replica.

    Notes
    -----
    On a CPU machine without ``boltz`` installed, this raises a
    ``RuntimeError`` with the SLURM template for the cluster.
    """

    binary = binary or _detect_boltz_binary()
    if binary is None:
        raise RuntimeError(
            "boltz CLI not found. This wrapper expects to run on a GPU node "
            "with `conda activate boltz` followed by `pip install boltz`. "
            "See docs/cluster.md for the SLURM template."
        )

    yaml_path = build_yaml_input(inputs, work_dir)
    out_dir = Path(work_dir) / f"{inputs.target_id}_{model}"
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        binary,
        "predict",
        str(yaml_path),
        "--out_dir",
        str(out_dir),
        "--diffusion_samples",
        str(num_samples),
        "--seed",
        str(inputs.seed),
    ]
    if model == "boltz2":
        cmd += ["--model", "boltz2"]
    if use_msa_server:
        cmd += ["--use_msa_server"]
    if extra_args:
        cmd += list(extra_args)

    log.info("$ %s", " ".join(cmd))
    subprocess.run(cmd, check=True, env={**os.environ})

    # Collect produced PDB / cif files from {out_dir}/{target_id}/
    target_subdir = out_dir / inputs.target_id
    samples = []
    for i, p in enumerate(sorted(target_subdir.glob("*.cif"))):
        samples.append(
            EnsembleSample(
                target_id=inputs.target_id,
                sample_id=i,
                structure_path=p,
            )
        )
    return TeacherOutputs(
        target_id=inputs.target_id,
        teacher_name=model,
        samples=samples,
        output_dir=out_dir,
    )


# ---- SLURM submission helper ----

SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=boltz-{target_id}
#SBATCH --gres=gpu:a100:1
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/boltz_{target_id}_%j.out

set -euo pipefail
source ~/.bashrc
conda activate boltz

uv run python -m boltzrna_diff.teachers.boltz_run \\
    --target {target_id} \\
    --rna {rna_sequence} \\
    --smiles "{ligand_smiles}" \\
    --model {model} \\
    --num-samples {num_samples} \\
    --out {out_dir}
"""


def render_slurm_script(
    target_id: str,
    rna_sequence: str,
    ligand_smiles: str,
    *,
    model: ModelName = "boltz1",
    num_samples: int = 16,
    out_dir: str = "runs/boltz",
) -> str:
    return SLURM_TEMPLATE.format(
        target_id=target_id,
        rna_sequence=rna_sequence,
        ligand_smiles=ligand_smiles,
        model=model,
        num_samples=num_samples,
        out_dir=out_dir,
    )
