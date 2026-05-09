"""Teacher-model wrappers for cofolding distillation.

Three GPU-only models are wrapped here:

* :mod:`boltzrna_diff.teachers.boltz` — Boltz-1 / Boltz-2 cofolding (primary
  teacher).
* :mod:`boltzrna_diff.teachers.rhofold` — RhoFold+ ensemble sampling for
  RNA-only conformations.
* :mod:`boltzrna_diff.teachers.af3` — AlphaFold 3 (or AF3-class) refinement
  on a small subset.

All actual inference happens on a GPU cluster; this VM is CPU-only. The
modules here therefore expose:

1. **Pure-Python adapters** that build inputs, validate them, and write them
   to disk in the format expected by each upstream tool.
2. **SLURM-friendly runner** functions that invoke the model binaries via
   ``subprocess`` and parse their outputs.
3. **Lazy imports** — heavy dependencies (boltz, openmm, etc.) are imported
   only inside :func:`run` so that the rest of the package can be installed
   with ``uv sync --extra dev`` on a CPU-only workstation.

See ``docs/cluster.md`` for the conda environment and SLURM submission
templates.
"""

from .types import EnsembleSample, TeacherInputs, TeacherOutputs

__all__ = ["EnsembleSample", "TeacherInputs", "TeacherOutputs"]
