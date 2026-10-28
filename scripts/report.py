"""
report.py
---------
Builds the design-rationale document from the pipeline's own outputs.

The report is generated from metrics.json rather than written by hand, so
the numbers in the text cannot drift away from the numbers the code
produced. Re-running the pipeline re-writes the report.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output"


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join("---" for _ in cols) + "|"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def build() -> Path:
    m = json.loads((OUT / "metrics.json").read_text())
    primers = pd.read_csv(OUT / "primers.csv")

    nai, opt = m["naive"], m["optimised"]
    fixture = m["vector_is_synthetic_fixture"]

    checks_df = pd.DataFrame([
        {"Check": c["name"],
         "Result": "PASS" if c["passed"] else "FAIL",
         "Detail": c["detail"]}
        for c in m["checks"]
    ])

    opt_df = pd.DataFrame([
        {"Metric": "Length (nt)", "Naive": nai["length_nt"], "Optimised": opt["length_nt"]},
        {"Metric": "CAI", "Naive": f"{nai['CAI']:.3f}", "Optimised": f"{opt['CAI']:.3f}"},
        {"Metric": "GC (%)", "Naive": f"{nai['GC_percent']:.1f}", "Optimised": f"{opt['GC_percent']:.1f}"},
        {"Metric": "Min 50 nt window GC (%)", "Naive": f"{nai['GC_window_min']:.0f}", "Optimised": f"{opt['GC_window_min']:.0f}"},
        {"Metric": "Max 50 nt window GC (%)", "Naive": f"{nai['GC_window_max']:.0f}", "Optimised": f"{opt['GC_window_max']:.0f}"},
        {"Metric": "Classic rare codons", "Naive": nai["rare_codons"], "Optimised": opt["rare_codons"]},
    ])

    primer_cols = ["Name", "Sequence (5'->3')", "Length", "Anneal Tm (C)",
                   "Anneal GC (%)", "Hairpin dG (kcal/mol)",
                   "Self-dimer dG (kcal/mol)", "Warnings"]
    primer_df = primers[[c for c in primer_cols if c in primers.columns]]

    warning_banner = ""
    if fixture:
        warning_banner = (
            "> **This run used a synthetic test fixture, not pET-28a(+).**\n"
            "> The assembly and verification code is exercised end to end, but the\n"
            "> vector is a stand-in. Re-run with `--vector` pointing at the real\n"
            "> GenBank file before treating any coordinate here as real.\n\n"
        )

    ep = m.get("expression_product", {})
    total_aa = ep.get("total_aa", "?")
    leader_aa = ep.get("leader_aa", "?")
    leader = ep.get("leader", "?")
    tail = ep.get("c_terminal_tail", "?")
    mass = ep.get("mass_kDa", "?")
    target_aa = m["protein_length_aa"]

    doc = f"""# Design of a codon-optimised RsTAL expression construct in pET-28a(+)

*In-silico design only. Nothing in this document has been built or tested at the bench.*

{warning_banner}## 1. What this is and why TAL

Tyrosine ammonia-lyase converts L-tyrosine into p-coumaric acid. That reaction is
the gateway into the phenylpropanoid pathway, which is where most plant
polyphenols, flavonoids and stilbenes come from. Any attempt to make those
compounds in a microbe starts with getting an ammonia-lyase to work in
*E. coli*, so TAL is a sensible first target rather than an arbitrary one.

The enzyme used here is the *Rhodobacter sphaeroides* enzyme (RsTAL,
UniProt Q3IWB0), chosen because it is structurally characterised: there are
crystal structures of the wild type bound to p-coumarate and of the H89F
variant, which gives a published reference point for later variant work.

**Target protein:** {m['protein_length_aa']} aa
**Provenance:** {m['protein_provenance']}

## 2. A sequence trap worth recording

The crystallographic chain is 521 residues, but one of those positions is not a
standard amino acid: it is MDO, the 4-methylideneimidazole-5-one (MIO)
prosthetic group, whose parent residues are Ala-Ser-Gly.

MIO is not added by the cell. It forms autocatalytically from the enzyme's own
backbone at that triad. So the gene has to encode Ala-Ser-Gly there, and the
coding sequence is {m['protein_length_aa']} residues, not 521. Taking the
structural chain at face value would mean ordering a gene with a two-residue
deletion in the middle of the catalytic centre, and the protein would express
and purify perfectly well while being completely dead.

## 3. Codon optimisation

Reverse translation and optimisation were done with DNA Chisel against the
Kazusa *E. coli* K-12 codon usage table, with forbidden restriction sites
treated as hard constraints in the same solve.

Doing it in one solve matters. Optimise first and strip sites afterwards and you
can reintroduce a rare codon; strip sites first and then optimise and you can
recreate the site.

{md_table(opt_df)}

Constraints applied:

- translation identical to the input protein (hard)
- no NdeI, XhoI, NcoI, NheI, BamHI, EcoRI, SacI, SalI, HindIII, NotI or EcoRV site
- GC within 30-70% in every 50 nt window, and 35-62% overall
- **5' ramp:** GC held to 25-50% across the first 15 codons
- **primer-binding windows:** no short hairpins in the first or last 60 nt; GC 45-62% across the last 45 nt

The 5' ramp is the constraint that is easiest to leave out and most expensive to
get wrong. A stable hairpin over the ribosome binding site and start codon blocks
translation initiation, and naively maximising CAI makes this *worse* for a
protein from a high-GC organism, because the preferred codons are themselves
GC-rich. The GC trace shows the deliberate dip at the 5' end.

Overall GC lands at {opt['GC_percent']:.0f}%, which is high. That is a property of
the protein, not a failure of the optimiser: RsTAL comes from a high-GC organism
and is rich in Ala, Arg, Gly and Pro, whose *E. coli*-preferred codons are
GC-rich. Pushing GC lower would mean giving up CAI.

![Codon metrics](codon_metrics.png)

![GC profile](gc_profile.png)

Forbidden sites remaining in the optimised insert: **{m['forbidden_sites_remaining'] or 'none'}**

## 4. Vector, cloning strategy, and what actually gets expressed

**Vector:** pET-28a(+), 5,369 bp. T7 promoter for high-level expression, *lac*
operator for IPTG induction, His6 for single-step IMAC, kanamycin resistance
(more stable in long cultures than ampicillin, where secreted beta-lactamase
degrades the drug and lets untransformed cells grow).

**Strategy:** NdeI / XhoI, stop codon omitted.

### Two things the real vector sequence settled

**The cassette is stored on the bottom strand.** In the SnapGene distribution of
pET-28a(+), the T7 promoter, start codon and MCS run in *decreasing* coordinates:
XhoI is at 159 and NdeI at 239. Digesting the stored strand as written would
insert the gene backwards, and every frame check would still pass, because the
arithmetic would be internally consistent. The pipeline detects the orientation
from the T7 promoter and reverse-complements before cutting.

**NdeI does not remove the N-terminal tag.** This is the opposite of what is
often assumed. The relevant stretch of pET-28a(+) reads:

```
AAGGAGATATACC ATG GGC AGC AGC CAT CAT CAT CAT CAT CAC AGC AGC GGC
              M   G   S   S   H   H   H   H   H   H   S   S   G
              CTG GTG CCG CGC GGC AGC CAT ATG
              L   V   P   R   G   S   H   M  <- the NdeI ATG
```

The NdeI site sits **downstream** of the His6 and thrombin site and **in frame**
with the vector's own start codon, which is the one the ribosome binding site
serves. So the product carries the N-terminal leader as well as the C-terminal
tag. Cloning into **NcoI** is what removes the N-terminal tag, and that is not
available here without mutating residue 2, since NcoI (CCATGG) requires the
second codon to begin with G and RsTAL residue 2 is Leu.

### The expressed product

| Property | Value |
|---|---|
| Length | {total_aa} aa |
| Composition | {leader_aa} aa leader + {target_aa} aa RsTAL + LEHHHHHH |
| N-terminal leader | `{leader}` |
| C-terminal tail | `{tail}` |
| Predicted mass | {mass} kDa (untagged RsTAL is 54.9 kDa) |

The leader is a mixed blessing. Earlier in the design I wanted to keep bulk away
from the N-terminus, because the MIO group forms within the protein's own chain
and an N-terminal extension sits closer to it than a C-terminal one. Twenty extra
residues is not nothing.

What makes it acceptable is the thrombin site in the middle of the leader.
Cleaving after LVPR removes the His6 and leaves GSHM ahead of the native
sequence, so the N-terminal bulk is removable after purification rather than
permanent. The two tags also give two independent purification handles, which is
useful when you do not yet know whether the protein behaves.

If expression or activity disappoints, the comparison to run is this construct
against an NcoI/XhoI version carrying an L2V substitution, which trades a
conservative point mutation for a clean N-terminus.

## 5. Primers

{md_table(primer_df)}

Two notes on how these numbers were produced.

Tm is the nearest-neighbour calculation from primer3, not the
4*(G+C) + 2*(A+T) rule, which is badly wrong above about 20 nt. For a primer
carrying a restriction tail, the Tm that governs the first PCR cycle is the Tm
of the **annealing portion only**, because the tail has nothing to anneal to
yet. That is the figure tabulated above, and the one to set the annealing
temperature from.

Secondary structure is judged on free energy, not melting temperature. primer3
reports a hairpin Tm that can look alarming for a 3 bp stem whose dG is only
about -2 kcal/mol, a structure that never forms at annealing temperature.
Thresholds here are dG < -3 kcal/mol for hairpins and dG < -6 kcal/mol for
dimers.

Cloning pair heterodimer dG: **{m['cloning_pair_heterodimer_dG']} kcal/mol**.

Each cloning primer carries a 6 base clamp 5' of the restriction site, because
most enzymes cut poorly at the very end of a fragment. The clamp sequence was
chosen by searching candidates and scoring them, not fixed in advance.

## 6. In-silico verification

Simulated build: amplicon {m['amplicon_bp']} bp, digested insert
{m['insert_fragment_bp']} bp, final plasmid {m['plasmid_bp']} bp.

{md_table(checks_df)}

**{sum(1 for c in m['checks'] if c['passed'])} of {len(m['checks'])} checks passed.**

The NdeI check is the one that catches real errors. NdeI leaves a TA overhang,
so the digested insert begins `TATG` and the backbone ends `CA`: the start codon
only exists once the two are ligated. The check confirms CATATG is actually
reconstituted rather than assuming it.

## 7. Confirmation plan for the bench

| Step | Method | Expectation |
|---|---|---|
| Insert present | Colony PCR with the cloning primers | Product at ~{m['amplicon_bp']} bp; empty vector gives no product or a much smaller one |
| Orientation and size | Diagnostic digest, NdeI + XhoI | Two bands: insert ~{m['insert_fragment_bp']} bp and backbone |
| Sequence | Sanger from T7 promoter and T7 terminator primers, plus the internal primers above | Full coverage; reads overlap |
| Expression | BL21(DE3), IPTG induction, SDS-PAGE | Band near the predicted mass |
| Activity | p-coumarate absorbance near 310 nm | Rate above the empty-vector control |

Sanger reads give roughly 700-900 clean bases, which is why internal primers are
needed: the T7 primers alone cannot cover a {m['optimised']['length_nt']} bp insert.

## 8. Limitations

- Design only. Nothing here has been synthesised, cloned, expressed or assayed.
- Codon optimisation predicts expression; it does not guarantee it. CAI correlates with expression across genes, but is a weak predictor for any single gene.
- Solubility was not modelled. Ammonia-lyases of this family are large and homotetrameric, and inclusion bodies are a realistic outcome.
- The hairpin and dimer calculations are two-state thermodynamic models, not a simulation of what happens in a PCR tube.
- The 5' ramp constraint is based on an established principle, but the specific window used here (15 codons, 25-50% GC) is a reasonable default, not a value optimised for this gene.
- No experimental validation of the MIO triad assignment beyond the structural annotation.

## 9. Reproducing this

```bash
pip install biopython python_codon_tables dnachisel primer3-py pandas matplotlib
python3 scripts/run_all.py --protein data/Q3IWB0.fasta --vector data/pET-28a.gb --seed {m['seed']}
```

DNA Chisel's solver is stochastic, so the run is seeded. Seed **{m['seed']}**
reproduces every sequence, primer and coordinate in this document exactly.

## References

- Watts KT, Mijts BN, Lee PC, Manning AJ, Schmidt-Dannert C. Discovery of a substrate selectivity switch in tyrosine ammonia-lyase. *Chem Biol* 2006;13:1317-26. PMID 17185227
- Louie GV, Bowman ME, Moffitt MC, Baiga TJ, Moore BS, Noel JP. Structural determinants and modulation of substrate specificity in phenylalanine-tyrosine ammonia-lyases. *Chem Biol* 2006;13:1327-38. PMID 17185228. PDB 2O6Y, 2O78, 2O7B, 2O7D, 2O7E, 2O7F
- Sharp PM, Li WH. The codon adaptation index. *Nucleic Acids Res* 1987;15:1281-95
- Zulkower V, Rosser S. DNA Chisel, a versatile sequence optimizer. *Bioinformatics* 2020
- Untergasser A et al. Primer3, new capabilities and interfaces. *Nucleic Acids Res* 2012;40:e115
"""

    path = OUT / "DESIGN_REPORT.md"
    path.write_text(doc)
    return path


if __name__ == "__main__":
    print(f"wrote {build()}")
