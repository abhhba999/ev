# Speccard: pre-miR-21

> Default RNA target #1 — best-characterized pre-miRNA for small molecule binding and the most precedented RIBOTAC substrate.

## Biology
- pre-miR-21 (HSA-MIR-21) is processed from a long primary transcript by Drosha (DGCR8) into a ~72-nt hairpin pre-miRNA, then exported and cleaved by Dicer to give mature miR-21-5p.
- miR-21 is one of the most upregulated oncomiRs across solid tumors (NSCLC, glioblastoma, hepatocellular carcinoma, pancreatic, breast, colorectal). Knockdown reverses chemoresistance and reduces tumor growth in many xenograft models.
- The Drosha cleavage site within pre-miR-21 contains a bulge motif that is the prototype small-molecule binding site (Disney lab targeting paradigm).

## Why pre-miR-21 is the right pick for our project
- **Most published rSM data**: at least 5 distinct chemotypes (Disney lab series, dimers, Targapremir-21, RIBOTAC version, dual-action RNA degrader). DC50 / IC50 values are public.
- **Tractable wet-lab readouts**: in-vitro Drosha processing assay (gel- or qPCR-based), miR-21 northern, downstream PDCD4/PTEN protein rescue, in-cell tumor proliferation (e.g. MDA-MB-231).
- **First RIBOTAC validated in vivo**: Costales et al. 2020 (PNAS, https://doi.org/10.1073/pnas.2003136117) demonstrated targeted degradation of pre-miR-21 in cells and mouse xenografts. This is our retrospective ground truth.

## Structures available
| Source | Accession | Note |
|---|---|---|
| HARIBOSS | none (as of 2025-04 dump) | pre-miR-21 hairpin binders are mostly characterized in solution; no co-crystal in current HARIBOSS |
| PDB | **6OF1**, **6OF2**, **6CK5** (related) | small-molecule + pre-miR hairpin family; we will pull related entries via SMARTBind workflow |
| Predicted | RhoFold+ + AF3 | we will generate K=8-32 conformers |
| Secondary structure | mfold + ViennaRNA | hairpin with internal loop / bulge near Drosha site |

## Reference compounds (for retrospective benchmark)
| Compound | Chemotype | Activity (target) | Reference |
|---|---|---|---|
| Targapremir-21 | benzimidazole + dimethylamine | ~200 nM in vitro Drosha block | Velagapudi et al., *Nat Chem Biol* 2014; doi:10.1038/nchembio.1452 |
| Targapremir-21 dimer | bivalent | low-nM, improved cell potency | Costales et al., *J Am Chem Soc* 2017 |
| RIBOTAC-pre-miR-21 | warhead + RNase L recruiter | first-gen RIBOTAC, in-cell + in-vivo activity | Costales et al., *PNAS* 2020; doi:10.1073/pnas.2003136117 |
| Dual-acting degrader | catalytic-cleavage hybrid | further potency gain | Yang et al., *J Am Chem Soc* 2024 |

## Wet-lab plan (with our collaborator)
1. **Retrospective**: regenerate the published binders in our generative model and confirm Vina/Boltz-2 affinity ranks them in the top 5%.
2. **Prospective (5-20 compounds)**: order generated compounds (ChemDiv / Enamine custom) → in-vitro Drosha processing assay → northern miR-21 in MDA-MB-231 → PDCD4 rescue Western.
3. **Stretch RIBOTAC**: 1-2 RIBOTAC conjugates with our binder + published Kethoxal-style RNase L recruiter (Costales linker).

## Open questions for collaborator
- Preferred cell line(s)? (MDA-MB-231 is canonical; the collaborator may have HCC827, A549, or PANC-1 lines on hand.)
- Is the SHAPE-MaP / in-cell SHAPE infrastructure available? (would let us sample the conformation ensemble experimentally.)
