# Speccard: HIV-1 TAR

> Default RNA target #3 — historically the longest-pursued RNA drug target, lots of co-crystals, well-defined pocket; an established benchmark for RNA SBDD and a clinically meaningful viral target.

## Biology
- TAR (Trans-Activation Response element) is a 59-nt stem-loop at the 5' of all HIV-1 transcripts. It binds the viral protein Tat + host P-TEFb complex to recruit elongation factors.
- Disrupting Tat-TAR shuts down HIV-1 transcription elongation. >40 years of medicinal chemistry; no approved drug yet but multiple advanced tool compounds.
- In our project, TAR is the "Goliath" benchmark target: structurally well-characterized, biologically meaningful, but historically very hard to drug selectively.

## Why HIV-1 TAR is the right pick
- **Lots of structures**: HIV-1 / HIV-2 TAR co-crystals in HARIBOSS (1AJU, 1AKX, 1ARJ, 1UTS, 1UUD, 1UUI, 2L8H, 6CMN, ...). Multiple chemotypes — argininamide, aminoglycosides, peptide mimetics, drug-like small molecules.
- **Cofolding stress test**: TAR has a UCU bulge with multiple metastable conformations. A perfect target to demonstrate the value of *ensemble-aware* generation vs single-structure baselines.
- **Selective benchmark**: assays like FRET (Tat-TAR displacement), TAR cassette reporter, and HIV-1 replication in cell are well-established.

## Structures available (from our HARIBOSS dump)
| PDB | Length | Method | Notes |
|---|---|---|---|
| 1AJU / 1AKX | 28-30 nt | NMR | argininamide complex |
| 1ARJ | 28 nt | NMR | argininamide |
| 1UTS / 1UUD / 1UUI | 28 nt | NMR | drug-like small molecules |
| 2L8H | — | NMR | newer chemotype |
| 6CMN / 6CSL | — | X-ray | drug-like fragments |

## Reference compounds (for retrospective benchmark)
| Compound | Activity | Reference |
|---|---|---|
| Argininamide | mM IC50 (canonical reference) | Aboul-ela et al., *J Mol Biol* 1995 |
| Neomycin / paromomycin | sub-μM | aminoglycoside class |
| TAR-Tat displacement small molecules | μM | Lind et al., Davis & Hamy series |
| WPL-43 / WPL-67 | nM | newer chemotypes; Olenyuk lab |

## Wet-lab plan
1. **Retrospective**: train on 80% of TAR co-crystals, hold out 20% scaffold-distinct chemotypes, ask the model to recover the held-out actives.
2. **Prospective**: 10-20 compounds → Tat-TAR FRET displacement → HIV-1 cassette luciferase reporter → cytotoxicity counter-screen.
3. **Stretch**: ensemble vs single-conformation ablation. We expect ensemble-aware model to win on TAR by a large margin because of its known UCU-bulge flexibility.
