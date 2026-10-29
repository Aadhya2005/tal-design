# RsTAL / pET-28a(+) in-silico construct design

A reproducible pipeline that designs an *E. coli* expression construct for
tyrosine ammonia-lyase from *Rhodobacter sphaeroides* (RsTAL, UniProt Q3IWB0):
codon optimisation, restriction-site clearance, primer design, simulated
digest and ligation, and automated verification of the finished plasmid.

**This is a design exercise. Nothing here has been built or tested at the bench.**

## Why TAL

TAL converts L-tyrosine into p-coumaric acid, the entry point into the
phenylpropanoid pathway that produces most plant polyphenols, flavonoids and
stilbenes. Making those compounds in a microbe starts with getting an
ammonia-lyase working in *E. coli*. RsTAL specifically, because it is
structurally characterised, including the H89F variant that flips substrate
preference from tyrosine to phenylalanine, which gives a published reference
point for later variant work.

## What it does

| Stage | Tool | Output |
|---|---|---|
| Codon optimisation + site clearance | DNA Chisel, Kazusa *E. coli* K-12 table | optimised CDS, CAI/GC metrics |
| Primer design | primer3 (nearest-neighbour thermodynamics) | cloning + sequencing primers, scored |
| Digest and ligation | Biopython `Bio.Restriction` | finished plasmid as GenBank |
| Verification | custom checks | 8-point pass/fail report |
| Independent audit | `audit.py`, re-reads the GenBank | product, tags, markers, digest |
| Figures and report | matplotlib | plasmid map, GC profile, design document |

## Running it

Requires Python 3.10 or newer. See `SETUP.md` for a full walkthrough,
including opening this folder in VS Code (a `.vscode/` config is included:
Run and Debug has `run_all.py` and `audit.py` pre-configured).

```bash
pip install -r requirements.txt

python3 scripts/run_all.py \
  --protein data/Q3IWB0.txt \
  --vector  data/pET-28a.dna \
  --seed    1

python3 scripts/report.py
python3 scripts/audit.py
```

File names above match `SETUP.md`. `run_all.py` detects format from the
extension, so `.txt`/`.dat`/`.sp` (Swiss-Prot) or `.fasta` both work for the
protein, and `.dna` (SnapGene), `.gb`/`.gbk`/`.genbank`, or `.fasta` all work
for the vector — use whatever your downloaded files are actually named.

`audit.py` deliberately does not reuse the build code's own checks. It re-reads
the written GenBank from disk and asks different questions, so a mistake in the
build logic cannot pass by being consistent with itself. It is what caught the
N-terminal tag error above.

Inputs are not committed. Fetch them yourself:

- **Protein** — UniProt `Q3IWB0`, download as Text (Swiss-Prot flat file) or FASTA
- **Vector** — pET-28a(+), from SnapGene's plasmid library (`.dna`) or Addgene (GenBank)

Omit `--vector` and the pipeline runs against a synthetic test fixture instead.
Every output from a fixture run is labelled as such. The fixture exists so the
assembly and verification code can be proven correct independently of the real
vector; it is not pET-28a(+) and its coordinates mean nothing.

## Five things the pipeline gets right that are easy to get wrong

**The MIO triad.** The RsTAL crystallographic chain is 521 residues, but one
position is MDO, the 4-methylideneimidazole-5-one cofactor, whose parents are
Ala-Ser-Gly. MIO forms autocatalytically from the protein's own backbone, so
the gene must encode all three residues and the real CDS is 523 codons. Reading
the structural chain literally gives a two-residue deletion in the catalytic
centre, and the resulting protein would express and purify normally while being
completely inactive.

**The 5' ramp.** The first ~15 codons of a bacterial CDS should be AT-biased and
structure-poor, because a hairpin over the ribosome binding site blocks
initiation. Naively maximising CAI makes this worse for a protein from a high-GC
organism, since the preferred codons are themselves GC-rich. The optimiser
applies an explicit GC constraint over that window; the dip is visible in
`gc_profile.png`.

**Vector strand orientation.** SnapGene's pET-28a(+) stores the expression
cassette on the bottom strand: the T7 promoter, start codon and MCS run in
decreasing coordinates, and XhoI (159) comes before NdeI (239). Digest the
stored strand as written and the insert goes in backwards, while every frame
check still passes because the arithmetic stays internally consistent. The
orientation is detected from the T7 promoter rather than assumed.

**NdeI does not strip the N-terminal tag.** Contrary to a common assumption,
the NdeI site in pET-28a(+) sits downstream of the His6 and thrombin site and
in frame with the vector's own start codon, which is the one the RBS serves.
The product is therefore `MGSSHHHHHHSSGLVPRGSH` + RsTAL + `LEHHHHHH`, 551 aa,
58.1 kDa, tagged at both ends. NcoI is the site that removes the N-terminal
tag, and it is unavailable here without mutating residue 2. The pipeline finds
the real translational start from the RBS instead of assuming the insert's own
ATG is used.

**dG, not Tm, for secondary structure.** primer3 reports a hairpin melting
temperature that can read above 60 °C for a 3 bp stem whose dG is only about
-2 kcal/mol, a structure that never forms at annealing temperature. Gating
primer warnings on Tm rejects perfectly good primers. The thresholds here are
dG < -3 kcal/mol for hairpins and dG < -6 kcal/mol for dimers.

## Reproducibility

DNA Chisel's solver is stochastic. Without a fixed seed the same input produces
a different valid sequence on every run, and the primers, coordinates and
reported CAI all shift with it. `--seed` is threaded through to numpy, so a
given seed reproduces every number in the report exactly.

The report itself is generated from `metrics.json`, not written by hand, so the
prose cannot drift away from what the code actually computed.


## Two fixes found by reusing this pipeline

Running these modules on a different protein (plant O-methyltransferases, in a
later project) turned up two faults that a single test case had hidden.

**A latent crash in the optimiser.** Two of the end-anchored AvoidHairpins
constraints trip a bug in DNA Chisel: localising a constraint returns None,
and the solver then dereferences it. Whether it fires depends on the sequence.
RsTAL never triggered it; three of the four 6OMTs did. optimise.py now tries
the full constraint set first and falls back to a reduced one only when that
fails, and it returns dropped_constraints, so a weakened run shows up instead
of passing silently.

**A primer check that flagged the wrong thing.** Both cloning sites are
palindromes, so any primer carrying CATATG or CTCGAG can pair with a copy of
itself through the site. On free energy alone every one of those primers looks
bad. But the dimers melt at around 5 C and nothing is left of them at a 60 C
anneal. check_primer now asks two questions, whether the dimer is
energetically strong and whether it is still stable near the annealing
temperature, and it logs the harmless case as a note instead of a warning.

Same lesson as the hairpin thresholds above: one number does not decide it. The
question is whether the structure holds together at the temperature the reaction
actually runs at.

Neither change alters this project’s output. The RsTAL run gives the same
sequence, the same primers and the same 8 of 8 checks as before.


## Layout

```
.vscode/               Run/Debug config and editor settings for VS Code
scripts/
  codon_tools.py       CAI, GC, rare-codon metrics
  optimise.py          DNA Chisel problem definition and site scanning
  primers.py           primer design, scoring and QC
  assemble.py          digest, ligation, verification checks
  figures.py           plasmid map, GC profile, metric bars
  make_test_vector.py  synthetic fixture (not pET-28a)
  run_all.py           driver
  report.py            builds DESIGN_REPORT.md from metrics.json
  audit.py             independent re-check of the finished construct
output/
  DESIGN_REPORT.md     the design-rationale document
  pET28a_RsTAL.gb      finished construct
  primers.csv          primer table with QC
  metrics.json         every computed number
  *.png                figures
```

## Limitations

Design only, no experimental validation. Codon optimisation correlates with
expression across many genes but predicts it weakly for any single one.
Solubility is not modelled, and ammonia-lyases of this family are large
homotetramers for which inclusion bodies are a realistic outcome. The
thermodynamic calculations are two-state models, not simulations of a PCR.

## References

- Watts KT et al. Discovery of a substrate selectivity switch in tyrosine ammonia-lyase. *Chem Biol* 2006;13:1317-26
- Louie GV et al. Structural determinants and modulation of substrate specificity in phenylalanine-tyrosine ammonia-lyases. *Chem Biol* 2006;13:1327-38
- Sharp PM, Li WH. The codon adaptation index. *Nucleic Acids Res* 1987;15:1281-95
- Zulkower V, Rosser S. DNA Chisel, a versatile sequence optimizer. *Bioinformatics* 2020
- Untergasser A et al. Primer3, new capabilities and interfaces. *Nucleic Acids Res* 2012;40:e115


## how i made this 


The code in this repository was written with AI assistance (Claude). I chose
the target enzyme, scoped the project, supplied the real UniProt and vector
sequence files, directed every design decision described below, and verified
the output myself — both by reading the generated plasmid in Benchling and by
running the included independent audit script. I did not hand-write the
optimisation, primer-design, or verification code.

I'm treating this as a starting point for learning the underlying molecular
biology and bioinformatics properly, not as a finished claim of hands-on lab
or coding experience.
