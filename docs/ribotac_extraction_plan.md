# RIBOTAC literature extraction plan

> Goal: ~200-400 RIBOTAC-class compounds with structured metadata, suitable as a leave-target-out benchmark and seed dataset for RIBOTAC-Diff.

## Schema

We follow `boltzrna_diff.data.schema.RibotacRecord`:

- `record_id`            — `ribotac_<author>_<year>_<n>`
- `source`               — `'ribotac_lit_v0'`
- `source_version`       — `'2026-04'`
- `citation`             — DOI of the paper
- `rna_target_name`      — e.g. `'pre-miR-21'`, `'pre-miR-96'`, `'JUN'`, `'MALAT1'`
- `rna_target_sequence`  — best-known sequence (Drosha-cleaved hairpin, full lncRNA segment, ...)
- `rna_class`            — RnaTargetClass enum
- `smiles`               — full RIBOTAC SMILES (drawn from SI or reproduced via ChemDraw → SMILES)
- `inchikey`             — auto-derived from SMILES
- `rna_binder_smiles`    — sub-SMILES of the RNA-binding warhead (when authors give it)
- `linker_smiles`        — when given
- `rnase_l_recruiter_smiles` — usually a known kethoxal- or 2'-5'A4-mimic moiety
- `recruiter_kind`       — string label
- `dc50_nM`              — half-maximal degradation concentration in nM
- `dmax_percent`         — max degradation observed
- `cell_line`            — e.g. `'MDA-MB-231'`, `'HCT-116'`, `'NSCLC'`
- `assay_kind`           — `'degradation_dc50'`

## Source list (priority order)

### Disney lab (the main RIBOTAC origin)
1. Costales et al., *PNAS* **2020** — first pre-miR-21 RIBOTAC, in vivo. doi:10.1073/pnas.2003136117
2. Liu et al., *J Am Chem Soc* **2021** — RIBOTAC for pre-miR-96.
3. Haniff et al., *Nat Chem Biol* **2020** — pre-miR-525-5p (HCC).
4. Tong et al., *Cell Chem Biol* **2022** — JUN mRNA.
5. Velagapudi et al., *J Am Chem Soc* **2022** — MALAT1 selective degradation.
6. Disney lab review, *Nat Rev Drug Discov* **2024** — comprehensive RIBOTAC retrospective.
7. Childs-Disney lab follow-up papers on pre-miR-21 (2022-2025) — chemotype expansion.
8. Yang et al., *J Am Chem Soc* **2024** — dual-acting RNA degraders.
9. Zhang et al., *Nat Commun* **2025** — TaRiboTAC tumor-microenvironment-activated RIBOTACs. doi:10.1038/s41467-025-56691-3

### Other sources to scan
- Arrakis Therapeutics published programs (DM1; ATX-1700-class).
- Expansion Therapeutics / Disney lab repeats programs.
- Skyhawk Therapeutics (mostly splicing modifiers, but bifunctional patents).
- Novartis / Roche / Pfizer patent literature 2022-2025 with `RIBOTAC` / `RNA degrader` keywords (USPTO + WO).
- RNA-degrader-related arXiv / bioRxiv 2024-2026 (search "RNase L recruiter", "RIBOTAC", "bifunctional RNA degrader").

## Workflow

For each source:
1. Pull PDF + SI.
2. Parse compound numbering vs full structure tables (often in SI).
3. Reproduce the structure in ChemDraw → SMILES (or use OSRA/Imago for OCR). Reviewer (lin) validates a 10% random sample.
4. For each compound: identify RNA target, cell line, DC50 / Dmax / IC50 / EC50 from main text Table or SI.
5. Append a row to `data/external/ribotac_v0.csv`.
6. Validate against `RibotacRecord` schema.

## Quality gates
- 100% of rows must round-trip SMILES through RDKit (`MolFromSmiles → MolToSmiles`).
- A second-pass LLM-assisted summary cross-checks DC50 unit + cell line from the article abstract.
- Each entry must have a working DOI.
- Ambiguous structures (e.g. patents with Markush) are excluded from v0; flagged for v1.

## Stretch (v1)
- Scrape WIPO / USPTO patents (RIBOTAC / RNA degrader) for the missing tail.
- Negative-control set: ~500 size- and chemotype-matched non-RIBOTAC bifunctionals (PROTAC, molecular glue) — keeps the model honest about what is or isn't an RNA degrader.
