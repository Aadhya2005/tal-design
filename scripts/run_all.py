"""
run_all.py
----------
Driver. Runs the whole design from inputs to deliverables.

    python3 scripts/run_all.py --protein data/Q3IWB0.fasta --vector data/pET-28a.gb

If --vector is omitted a synthetic test fixture is used instead, and every
output is clearly marked as a fixture run.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))

from assemble import (annotate, build_construct, digest_insert,
                      find_expression_product, pcr_product, verify)
from codon_tools import summarise, translate_cds
from figures import circular_map, codon_usage_bars, gc_profile
from optimise import find_sites, optimise_cds
from primers import (T7_PROMOTER_PRIMER, T7_TERMINATOR_PRIMER,
                     design_cloning_primers, design_internal_seq_primers,
                     pair_heterodimer_dg, primer_table_rows, Primer,
                     check_primer)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output"


def load_protein(path: Path | None) -> tuple[str, str]:
    """Return (sequence, provenance_note)."""
    if path is not None and path.suffix.lower() in (".txt", ".dat", ".sp"):
        rec = next(SeqIO.parse(str(path), "swiss"))
        return str(rec.seq).strip().upper(), (
            f"UniProtKB/Swiss-Prot {rec.id} ({rec.annotations.get('entry_name','')}), "
            f"sequence version {rec.annotations.get('sequence_version','?')}"
        )
    if path is None:
        prov = ROOT / "data" / "rstal_from_pdb_UNVERIFIED.txt"
        seq = prov.read_text().strip().replace("X", "ASG")
        return seq, (
            "PROVISIONAL: chain sequence read from the RCSB entry for the "
            "RsTAL structures, with the MDO modified residue expanded back to "
            "its Ala-Ser-Gly parents. Not verified against UniProt."
        )
    rec = next(SeqIO.parse(str(path), "fasta"))
    return str(rec.seq).strip().upper(), f"UniProt FASTA: {rec.id}"


T7_PROMOTER_CORE = "TAATACGACTCACTATAG"


def load_vector(path: Path | None):
    """
    Load the vector and make sure the expression cassette reads 5'->3' on the
    top strand.

    This is not cosmetic. SnapGene's pET-28a(+) file stores the cassette on the
    bottom strand: the T7 promoter, start codon and MCS run in *decreasing*
    coordinates, and XhoI (159) comes before NdeI (239). Digesting the stored
    strand as written would put the insert in backwards, downstream of nothing,
    and every frame check would still pass because the arithmetic is internally
    consistent. The orientation is therefore detected from the T7 promoter
    rather than assumed.
    """
    if path is None:
        from make_test_vector import build
        return build(), True, False

    suffix = path.suffix.lower()
    if suffix == ".dna":
        fmt = "snapgene"
    elif suffix in (".gb", ".gbk", ".genbank"):
        fmt = "genbank"
    else:
        fmt = "fasta"

    rec = next(SeqIO.parse(str(path), fmt))
    rec.annotations.setdefault("molecule_type", "DNA")

    top = str(rec.seq).upper()
    bottom = str(rec.seq.reverse_complement()).upper()
    flipped = False
    if T7_PROMOTER_CORE not in top:
        if T7_PROMOTER_CORE in bottom:
            rc = rec.reverse_complement()
            rc.id = rec.id if rec.id != "<unknown id>" else path.stem
            rc.name = path.stem
            rc.description = f"{path.stem} (reverse-complemented so the T7 cassette reads on the top strand)"
            rc.annotations["molecule_type"] = "DNA"
            rc.annotations["topology"] = "circular"
            rec = rc
            flipped = True
        else:
            raise ValueError(
                "T7 promoter not found on either strand; is this a pET vector?"
            )

    if rec.id in ("", "<unknown id>"):
        rec.id = path.stem
    rec.annotations.setdefault("topology", "circular")
    return rec, False, flipped


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", type=Path, default=None)
    ap.add_argument("--vector", type=Path, default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--name", default="pET28a_RsTAL")
    args = ap.parse_args()

    OUT.mkdir(exist_ok=True)

    protein, prov_note = load_protein(args.protein)
    vector, is_fixture, flipped = load_vector(args.vector)

    print(f"Protein      : {len(protein)} aa")
    print(f"Provenance   : {prov_note}")
    print(f"Vector       : {vector.id}, {len(vector.seq)} bp"
          + ("  [SYNTHETIC FIXTURE]" if is_fixture else ""))
    if flipped:
        print("               reverse-complemented: the T7 cassette was stored "
              "on the bottom strand")
    print(f"Residue 89   : {protein[88]} (expect H, the TAL specificity determinant)")

    # --- optimise ---------------------------------------------------------
    print("\nOptimising codons...")
    res = optimise_cds(protein, seed=args.seed)
    naive, opt = res["naive_cds"], res["optimised_cds"]
    assert translate_cds(opt) == protein, "translation changed during optimisation"

    naive_stats = summarise(naive, "naive")
    opt_stats = summarise(opt, "optimised")
    print(f"  CAI  {naive_stats['CAI']:.3f} -> {opt_stats['CAI']:.3f}")
    print(f"  GC   {naive_stats['GC_percent']:.1f}% -> {opt_stats['GC_percent']:.1f}%")
    print(f"  forbidden sites remaining: {find_sites(opt) or 'none'}")

    # --- primers ----------------------------------------------------------
    print("\nDesigning primers...")
    fwd, rev = design_cloning_primers(opt)
    internals = design_internal_seq_primers(opt)
    t7f = check_primer(Primer("T7_promoter", T7_PROMOTER_PRIMER, T7_PROMOTER_PRIMER,
                              "standard vector sequencing primer, reads into the 5' end",
                              reference=True))
    t7r = check_primer(Primer("T7_terminator", T7_TERMINATOR_PRIMER, T7_TERMINATOR_PRIMER,
                              "standard vector sequencing primer, reads into the 3' end",
                              reference=True))
    all_primers = [fwd, rev] + internals + [t7f, t7r]
    het = pair_heterodimer_dg(fwd, rev)
    print(f"  cloning pair heterodimer dG: {het:.1f} kcal/mol")
    for p in all_primers:
        flag = "OK" if not p.warnings else "; ".join(p.warnings)
        print(f"  {p.name:<16} {len(p.sequence):>2} nt  Tm {p.tm_anneal:5.1f}  {flag}")

    # --- assemble ---------------------------------------------------------
    print("\nAssembling...")
    amplicon = pcr_product(opt, fwd.sequence, rev.sequence)
    fragment = digest_insert(amplicon)
    construct = build_construct(vector, fragment, name=args.name)
    construct = annotate(construct, len(protein))
    print(f"  amplicon {len(amplicon)} bp -> insert {len(fragment)} bp "
          f"-> plasmid {len(construct.seq)} bp")

    # --- verify -----------------------------------------------------------
    print("\nVerification:")
    checks = verify(construct, protein, vector)
    for c in checks:
        print("  " + c.line())
    all_pass = all(c.passed for c in checks)
    print(f"\n  {sum(c.passed for c in checks)}/{len(checks)} checks passed")

    # --- write outputs ----------------------------------------------------
    (OUT / "RsTAL_optimised_CDS.fasta").write_text(
        f">RsTAL_codon_optimised_Ecoli len={len(opt)}nt seed={args.seed}\n"
        + "\n".join(opt[i:i + 60] for i in range(0, len(opt), 60)) + "\n"
    )
    SeqIO.write(construct, str(OUT / f"{args.name}.gb"), "genbank")

    rows = primer_table_rows(all_primers)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "primers.csv", index=False)

    metrics = {
        "protein_length_aa": len(protein),
        "protein_provenance": prov_note,
        "vector": vector.id,
        "vector_is_synthetic_fixture": is_fixture,
        "vector_reverse_complemented": flipped,
        "vector_length_bp": len(vector.seq),
        "seed": args.seed,
        "naive": {k: v for k, v in naive_stats.items() if k != "rare_detail"},
        "optimised": {k: v for k, v in opt_stats.items() if k != "rare_detail"},
        "forbidden_sites_remaining": find_sites(opt),
        "cloning_pair_heterodimer_dG": round(het, 2),
        "amplicon_bp": len(amplicon),
        "insert_fragment_bp": len(fragment),
        "plasmid_bp": len(construct.seq),
        "checks": [{"name": c.name, "passed": c.passed, "detail": c.detail}
                   for c in checks],
        "all_checks_passed": all_pass,
    }
    expr = find_expression_product(str(construct.seq))
    if expr.get("found"):
        from Bio.SeqUtils import molecular_weight
        from Bio.Seq import Seq as _Seq
        metrics["expression_product"] = {
            "leader": expr["leader"],
            "leader_aa": len(expr["leader"]),
            "total_aa": expr["length_aa"],
            "c_terminal_tail": expr["product"][-8:],
            "mass_kDa": round(molecular_weight(_Seq(expr["product"]),
                                               seq_type="protein") / 1000.0, 1),
            "start_bp_upstream_of_insert": expr["start_offset_before_insert"],
        }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2))

    # --- figures ----------------------------------------------------------
    print("\nFigures...")
    cds_len = len(protein) * 3
    plasmid_len = len(construct.seq)
    lead_bp = expr["start_offset_before_insert"] if expr.get("found") else 0
    feats = [
        ("RsTAL CDS (codon-optimised)", 1, cds_len, "#1f6f8b"),
        ("XhoI", cds_len + 1, cds_len + 6, "#c75b39"),
        ("C-terminal His6", cds_len + 7, cds_len + 24, "#4a7c59"),
        ("vector backbone", cds_len + 25, plasmid_len - lead_bp, "#b9b5ad"),
        ("N-term His6 + thrombin (retained)", plasmid_len - lead_bp + 1,
         plasmid_len, "#8a6fa8"),
    ]
    circular_map(plasmid_len, feats,
                 args.name + ("  [FIXTURE]" if is_fixture else ""),
                 str(OUT / "plasmid_map.png"))
    gc_profile(naive, opt, str(OUT / "gc_profile.png"))
    codon_usage_bars(naive_stats, opt_stats, str(OUT / "codon_metrics.png"))

    import report
    report.build()

    print(f"\nWrote outputs to {OUT}")
    for f in sorted(OUT.iterdir()):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
