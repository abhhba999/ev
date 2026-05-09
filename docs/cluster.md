# Cluster setup (GPU box, ≤8×A100)

This file lives outside `pyproject.toml` because Boltz-1, Boltz-2, AlphaFold 3, RhoFold+, OpenMM and openff-toolkit have native / conda-only dependencies that break uv resolution on a CPU-only workstation.

## 1. Base environment (mamba/conda)

```bash
mamba create -n boltzrna python=3.11 -y
mamba activate boltzrna

# CUDA-matched torch first (adjust cu version to your cluster).
mamba install -y pytorch=2.3 pytorch-cuda=12.1 -c pytorch -c nvidia

mamba install -y -c conda-forge openmm=8.1 mdtraj=1.10 openff-toolkit-base=0.16

# the boltzrna-diff project itself (pure Python deps via uv).
pip install -e '.[gpu]'
```

## 2. Boltz-1 / Boltz-2 (main teacher / refinement teacher)

```bash
git clone https://github.com/jwohlwend/boltz.git
cd boltz
pip install -e .

# weights (auto-downloaded on first run, ~25 GB):
boltz predict --help
```

## 3. RhoFold+ (RNA structure ensemble generation)

```bash
git clone https://github.com/ml4bio/RhoFold.git
cd RhoFold
pip install -e .
# weights:
wget https://example.org/rhofold_pretrained.pt -O pretrained/RhoFold_pretrained.pt
```

## 4. AlphaFold 3 (optional refinement teacher — academic use)

Follow https://github.com/google-deepmind/alphafold3 for license + weights.

## 5. SLURM template

```bash
sbatch scripts/slurm/train_student.sbatch
```

```sbatch
# scripts/slurm/train_student.sbatch
#SBATCH -J boltzrna-student
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:8
#SBATCH -c 64
#SBATCH --mem=512G
#SBATCH -t 7-00:00
#SBATCH -o logs/%x-%j.log

mamba activate boltzrna
python -m boltzrna_diff.training.run \
    +experiment=ensemble_distill_v1 \
    data=rna_rsm_bench_v0 \
    teacher=boltz1 \
    refiner=boltz2
```

## 6. Approx GPU budgets (Boltz-1 main teacher, Boltz-2 refinement)

| Step | Hardware | Time |
|---|---|---|
| RhoFold+ ensemble for 1 RNA | 1×A100 40 GB | 5-10 min |
| Boltz-1 cofolding 1 (RNA, ligand) | 1×A100 80 GB | 30-90 s |
| Boltz-2 cofolding 1 (RNA, ligand) | 1×A100 80 GB | 3-8 min |
| Student warm-up (CrossDocked2020) | 4×A100 | 3-5 d |
| RNA fine-tune + distillation | 4×A100 | 5-10 d |
| Full ensemble-aware run | 8×A100 | 7-14 d |

The whole training plan should fit inside ~30k GPU·h, comfortably under our budget. If we hit limits, the `teacher.cofolding_subsample` Hydra knob drops the per-batch teacher fraction from 100% to 25%.
