# BoltzRNA-Diff

> **Cofolding-distilled, conformation-aware generative design for RNA-targeted small molecules and RIBOTACs.**
>
> Target venue: ICLR 2027 (main track, deadline 2026-09-24)
> Plan-A' venue: *Nature Computational Science* / *Nature Methods*
> Default RNA targets (W1): pre-miR-21, FMN riboswitch, HIV-1 TAR

---

## What this repo will do

1. **RNA-rSM-Bench** — first unified benchmark for RNA-targeted small molecule design.
   Aggregates HARIBOSS, ROBIN, R-BIND 2.0, Inforna, RNAmigos2-supp, SMRTnet-supp into a
   common parquet schema with four OOD splits (scaffold, temporal, family, binding-site).
2. **Cofolding distillation** — uses Boltz-1 (main teacher) and Boltz-2 (refinement)
   as differentiable teachers to inject large-model priors into a small E(3)-equivariant
   flow-matching ligand generator that we can actually train on ≤8×A100.
3. **Conformation-ensemble conditioning** — RhoFold+ / AF3 + short MD ensembles per RNA
   target, plus an ensemble-attention generator that optimizes worst-case affinity over the
   conformation set rather than a single static structure.
4. **RIBOTAC-Bench / RIBOTAC-Diff** — first public benchmark and generative model for
   RNase L–recruiting RNA degraders (the "PROTAC for RNA" modality).

## Repo layout

```
boltzrna-diff/
├── pyproject.toml
├── README.md
├── configs/             # Hydra configs (data, model, split, training)
├── data/
│   ├── raw/             # untouched downloads
│   ├── processed/       # unified parquet schema
│   └── external/        # 3rd-party resources (HARIBOSS, RiboBind, ...)
├── docs/                # design notes, target speccards, methodology
├── notebooks/           # exploratory + figure generation
├── src/boltzrna_diff/
│   ├── data/            # loaders, schema, splits
│   ├── models/          # flow-matching student + cofolding-teacher wrappers
│   ├── training/        # train loops, distillation losses
│   ├── evaluation/      # metrics, retrospective hit recovery
│   └── scripts/         # CLI entrypoints (download_*, build_bench_*, run_*, ...)
└── tests/
```

## Hardware policy

Devin's workstation is CPU-only; this repo is structured so that:

- All **data acquisition, schema, splits, evaluation harness, RIBOTAC literature
  extraction, and unit tests** are pure Python / RDKit / pandas → run on any laptop.
- Anything that requires a GPU (Boltz-1/2 inference, RhoFold+, MD, training the flow
  matching student) is gated behind `[gpu]` / `[boltz]` extras and shipped as
  reproducible scripts plus a SLURM template (`docs/cluster.md`).

## Install

```bash
# CPU-only (data / schema / evaluation work)
uv pip install -e '.[dev]'

# Training / inference machine
uv pip install -e '.[gpu,boltz,md,dev]'
```

## Status (W1-W2)

- [x] Repo scaffold, packaging, lint config
- [ ] HARIBOSS loader + schema
- [ ] ROBIN, R-BIND, Inforna, RNAmigos2-supp, SMRTnet-supp loaders
- [ ] Four OOD split implementations
- [ ] RIBOTAC literature extraction v0
- [ ] Boltz-1 / Boltz-2 / RhoFold+ inference wrappers
- [ ] Speccards for default 3 RNA targets
