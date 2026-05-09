# Speccard: FMN riboswitch

> Default RNA target #2 — best-characterized riboswitch for ligand discovery; "PDBbind-style" benchmark target where structure is plentiful.

## Biology
- FMN (flavin mononucleotide) riboswitch is found in the 5'-UTR of bacterial flavin biosynthesis genes (e.g. *ribB*, *ribDEAH*).
- Native ligand FMN binds the aptamer with low-nM affinity and stabilizes a transcriptional terminator → represses gene expression.
- Antibacterial program: ribocil class (Merck) reached preclinical stage; the riboswitch is essential in *E. coli* and *B. subtilis*.

## Why FMN riboswitch is the right pick
- **Structural data is rich**: 10+ PDB co-crystal entries with diverse ligands in HARIBOSS (1FMN, 3F2T/W/X/Y, 3F4G/H, 5C45, 5KX9, 6BFB, ...). This is **the** well-mined RNA pocket for benchmarking.
- **Multiple chemotypes**: native FMN, reduced FMN, ribocil-A/B/C series, fluoroquinazoline scaffolds — gives us a clean *prospective scaffold-hop* test.
- **Wet-lab readouts are easy**: thermal melt (UV-Vis at 260 nm), 2-aminopurine fluorescence, ITC, *E. coli ribB-lacZ* reporter assay.

## Structures available (in HARIBOSS)
Confirmed from our `data/processed/hariboss.parquet`:
| PDB | Length | Method | Notes |
|---|---|---|---|
| 1FMN | 35 nt | NMR | aptamer core + FMN |
| 3F2T / 3F2W / 3F2X / 3F2Y | 107 nt | X-ray | full aptamer + FMN/reduced FMN |
| 3F4G / 3F4H | 110 nt | X-ray | aptamer + FMN/lumiflavin variants |
| 5C45 / 5KX9 | 109-110 nt | X-ray | ribocil-class hits |
| 6BFB | 110 nt | X-ray | fluoroquinazoline class |
| 6WJR | — | X-ray | newer chemotype |

## Reference compounds (for retrospective benchmark)
| Compound | Activity | Reference |
|---|---|---|
| FMN (native) | Kd ≈ 5 nM | Serganov et al., *Nature* 2009 |
| Ribocil-A/B/C | nM-μM in *E. coli* | Howe et al., *Nature* 2015; doi:10.1038/nature15542 |
| Ribocil derivatives | series of ~30 SAR points | Howe et al., *PNAS* 2016; doi:10.1073/pnas.1602393113 |

## Wet-lab plan
1. **Retrospective**: scaffold-hop test — train on 1FMN + 3F* class, hold out ribocil class, ask the model to rediscover ribocil-like chemotypes.
2. **Prospective**: 10-20 generated compounds → 2AP fluorescence + ITC → *ribB-lacZ* reporter for orthogonal in-cell readout.
3. **Stretch**: design selective binders that discriminate FMN riboswitch from RFN box homologs and from human FMN-binding proteins.
