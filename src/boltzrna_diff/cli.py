"""CLI entrypoints for boltzrna-diff."""

from __future__ import annotations

import logging
from pathlib import Path

import click

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")


@click.group()
def main() -> None:
    """boltzrna-diff command-line interface."""


@main.group("data")
def data_group() -> None:
    """Data acquisition and benchmark assembly."""


@data_group.command("download-hariboss")
@click.option(
    "--out",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("data/raw/hariboss"),
    show_default=True,
)
def download_hariboss(out: Path) -> None:
    """Fetch the HARIBOSS compounds + complexes CSVs."""

    from boltzrna_diff.data.hariboss import download

    a, b = download(out)
    click.echo(f"compounds → {a}\ncomplexes → {b}")


@data_group.command("build-hariboss")
@click.option(
    "--raw",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("data/raw/hariboss"),
    show_default=True,
)
@click.option(
    "--out",
    type=click.Path(dir_okay=False, path_type=Path),
    default=Path("data/processed/hariboss.parquet"),
    show_default=True,
)
def build_hariboss(raw: Path, out: Path) -> None:
    """Parse HARIBOSS CSVs into a unified parquet with the StructureRecord schema."""

    from boltzrna_diff.data.hariboss import load, to_records

    df = load(raw)
    n_valid = sum(1 for _ in to_records(df))
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out)
    click.echo(f"wrote {out} | rows={len(df)} | valid={n_valid}")


@data_group.command("download-robin")
@click.option(
    "--out",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("data/raw/robin"),
    show_default=True,
)
def download_robin(out: Path) -> None:
    """Sparse-clone the ROBIN GitHub repo and copy SMM_full_results CSVs."""

    from boltzrna_diff.data.robin import download

    p = download(out)
    click.echo(f"ROBIN CSVs in {p}")


@data_group.command("build-robin")
@click.option(
    "--raw",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("data/raw/robin"),
    show_default=True,
)
@click.option(
    "--out-binary",
    type=click.Path(dir_okay=False, path_type=Path),
    default=Path("data/processed/robin_binary.parquet"),
    show_default=True,
)
@click.option(
    "--out-per-target",
    type=click.Path(dir_okay=False, path_type=Path),
    default=Path("data/processed/robin_per_target.parquet"),
    show_default=True,
)
def build_robin(raw: Path, out_binary: Path, out_per_target: Path) -> None:
    """Parse ROBIN into BinaryBinderRecord + per-target InteractionRecord parquets."""

    from boltzrna_diff.data.robin import load_binary, load_per_target

    bdf = load_binary(raw)
    out_binary.parent.mkdir(parents=True, exist_ok=True)
    bdf.to_parquet(out_binary)
    click.echo(f"binary → {out_binary} | rows={len(bdf)} | rna_binders={int(bdf['is_rna_binder'].sum())}")

    pdf = load_per_target(raw)
    pdf.to_parquet(out_per_target)
    click.echo(
        f"per-target → {out_per_target} | rows={len(pdf)} | targets={pdf['rna_sequence'].nunique()} | actives={int(pdf['is_active'].sum())}"
    )


@data_group.command("describe")
@click.argument("parquet", type=click.Path(exists=True, dir_okay=False, path_type=Path))
def describe(parquet: Path) -> None:
    """Print a quick distribution summary for any RNA-rSM-Bench parquet."""

    import pandas as pd

    df = pd.read_parquet(parquet)
    click.echo(f"rows: {len(df)}")
    if "rna_class" in df.columns:
        click.echo("\nrna_class:")
        click.echo(df["rna_class"].value_counts().to_string())
    if "structural_method" in df.columns:
        click.echo("\nstructural_method:")
        click.echo(df["structural_method"].value_counts().to_string())
    if "ligand_inchikey" in df.columns:
        click.echo(f"\nunique ligands: {df['ligand_inchikey'].nunique()}")
    if "pdb_id" in df.columns:
        click.echo(f"unique PDB ids: {df['pdb_id'].nunique()}")


if __name__ == "__main__":
    main()
