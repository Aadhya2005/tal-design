"""
codon_tools.py
--------------
Sequence metrics for the RsTAL / pET-28a(+) construct design.

Everything here is computed from a published E. coli codon usage table
(Kazusa, E. coli K-12 W3110, taxid 316407) supplied by the
`python_codon_tables` package. No hand-entered frequency values.

CAI is implemented as Sharp & Li (1987):
    w(codon) = f(codon) / f(most frequent synonymous codon)
    CAI      = geometric mean of w over all codons in the CDS
Met and Trp have a single codon each (w = 1 always) and are excluded,
as are stop codons, which is the standard convention.
"""

from __future__ import annotations

import math
from collections import defaultdict

import python_codon_tables as pct
from Bio.Data import CodonTable
from Bio.Seq import Seq

ECOLI_TABLE_NAME = "e_coli_316407"

# Codons with low usage in E. coli that are conventionally flagged.
# These are reported, not auto-defined by a threshold, so the list is explicit.
CLASSIC_RARE_CODONS = {
    "AGG": "Arg", "AGA": "Arg", "CGA": "Arg", "CGG": "Arg",
    "CTA": "Leu",
    "ATA": "Ile",
    "CCC": "Pro",
    "GGA": "Gly",
}

_STANDARD = CodonTable.unambiguous_dna_by_id[1]


def load_ecoli_usage() -> dict[str, dict[str, float]]:
    """Return {aa: {codon: relative_frequency_within_that_aa}}."""
    return pct.get_codons_table(ECOLI_TABLE_NAME)


def relative_adaptiveness(usage: dict[str, dict[str, float]]) -> dict[str, float]:
    """Sharp & Li w values: codon frequency divided by the best synonymous codon."""
    w: dict[str, float] = {}
    for aa, codons in usage.items():
        best = max(codons.values())
        if best <= 0:
            continue
        for codon, freq in codons.items():
            w[codon] = freq / best
    return w


def codons_of(cds: str) -> list[str]:
    cds = cds.upper().replace("U", "T")
    if len(cds) % 3:
        raise ValueError(f"CDS length {len(cds)} is not a multiple of 3")
    return [cds[i:i + 3] for i in range(0, len(cds), 3)]


def cai(cds: str, w: dict[str, float] | None = None) -> float:
    """Codon Adaptation Index. Excludes Met, Trp and stop codons."""
    if w is None:
        w = relative_adaptiveness(load_ecoli_usage())

    single_codon_aas = {"M", "W"}
    logs: list[float] = []
    for codon in codons_of(cds):
        if codon in _STANDARD.stop_codons:
            continue
        aa = _STANDARD.forward_table.get(codon)
        if aa is None or aa in single_codon_aas:
            continue
        value = w.get(codon)
        if value is None or value <= 0:
            continue
        logs.append(math.log(value))

    if not logs:
        return float("nan")
    return math.exp(sum(logs) / len(logs))


def gc_percent(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return float("nan")
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def gc_windows(seq: str, window: int = 50, step: int = 10) -> list[tuple[int, float]]:
    """Sliding-window GC. Returns [(1-based window start, GC%)]."""
    seq = seq.upper()
    out: list[tuple[int, float]] = []
    if len(seq) < window:
        return [(1, gc_percent(seq))]
    for start in range(0, len(seq) - window + 1, step):
        out.append((start + 1, gc_percent(seq[start:start + window])))
    return out


def gc_extremes(seq: str, window: int = 50, step: int = 10) -> dict[str, float]:
    vals = gc_windows(seq, window, step)
    gcs = [g for _, g in vals]
    lo_pos, lo = min(vals, key=lambda t: t[1])
    hi_pos, hi = max(vals, key=lambda t: t[1])
    return {
        "window": window,
        "min_gc": lo, "min_at": lo_pos,
        "max_gc": hi, "max_at": hi_pos,
        "mean_gc": sum(gcs) / len(gcs),
    }


def rare_codon_report(cds: str) -> dict:
    """Count classic rare codons and report their 1-based codon positions."""
    found: dict[str, list[int]] = defaultdict(list)
    all_codons = codons_of(cds)
    for i, codon in enumerate(all_codons, start=1):
        if codon in CLASSIC_RARE_CODONS:
            found[codon].append(i)
    total = sum(len(v) for v in found.values())
    return {
        "total_rare": total,
        "n_codons": len(all_codons),
        "percent_rare": 100.0 * total / len(all_codons) if all_codons else float("nan"),
        "by_codon": {k: v for k, v in sorted(found.items())},
    }


def translate_cds(cds: str, to_stop: bool = False) -> str:
    return str(Seq(cds).translate(table=1, to_stop=to_stop))


def summarise(cds: str, label: str, w: dict[str, float] | None = None) -> dict:
    if w is None:
        w = relative_adaptiveness(load_ecoli_usage())
    rare = rare_codon_report(cds)
    ext = gc_extremes(cds)
    return {
        "label": label,
        "length_nt": len(cds),
        "length_codons": len(cds) // 3,
        "CAI": cai(cds, w),
        "GC_percent": gc_percent(cds),
        "GC_window_min": ext["min_gc"],
        "GC_window_max": ext["max_gc"],
        "GC_window_max_at": ext["max_at"],
        "rare_codons": rare["total_rare"],
        "rare_percent": rare["percent_rare"],
        "rare_detail": rare["by_codon"],
    }


if __name__ == "__main__":
    usage = load_ecoli_usage()
    w = relative_adaptiveness(usage)
    print(f"Loaded codon table: {ECOLI_TABLE_NAME}")
    print(f"Amino acid groups : {len(usage)}")
    print(f"Codons with w     : {len(w)}")
    print("\nSanity check - w must equal 1.0 for the top codon of each amino acid:")
    for aa in sorted(usage)[:6]:
        best = max(usage[aa], key=usage[aa].get)
        print(f"  {aa}: best codon {best}  w={w[best]:.3f}  freq={usage[aa][best]:.3f}")
