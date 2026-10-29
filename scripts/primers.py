"""
primers.py
----------
Design and quality-check the cloning and sequencing primers.

Tm is the nearest-neighbour calculation from primer3 (SantaLucia
thermodynamics), not the 4*(G+C) + 2*(A+T) rule of thumb, which is badly
wrong above about 20 nt.

For a cloning primer carrying a restriction tail, the Tm that governs the
FIRST PCR cycle is the Tm of the annealing portion only, because the tail
has nothing to anneal to yet. The full-length Tm applies from cycle three
onward. Both are reported; the annealing Tm is the one to set the
annealing temperature from.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import primer3

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

# Standard pET sequencing primers (Novagen). These read into the insert
# from upstream and downstream of the MCS.
T7_PROMOTER_PRIMER = "TAATACGACTCACTATAGGG"
T7_TERMINATOR_PRIMER = "GCTAGTTATTGCTCAGCGG"


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def gc_percent(seq: str) -> float:
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else float("nan")


@dataclass
class Primer:
    name: str
    sequence: str
    anneal_part: str
    purpose: str
    tail_note: str = ""
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    # Standard vector primers are published fixed sequences used under
    # sequencing conditions, not PCR. Holding them to PCR Tm windows would
    # report a failure that is not one.
    reference: bool = False

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def tm_anneal(self) -> float:
        return primer3.calc_tm(self.anneal_part)

    @property
    def tm_full(self) -> float:
        return primer3.calc_tm(self.sequence)

    @property
    def gc(self) -> float:
        return gc_percent(self.anneal_part)

    @property
    def hairpin_tm(self) -> float:
        return primer3.calc_hairpin(self.sequence).tm

    @property
    def hairpin_dg(self) -> float:
        """Free energy in kcal/mol (primer3 reports cal/mol)."""
        return primer3.calc_hairpin(self.sequence).dg / 1000.0

    @property
    def homodimer_tm(self) -> float:
        return primer3.calc_homodimer(self.sequence).tm

    @property
    def homodimer_dg(self) -> float:
        return primer3.calc_homodimer(self.sequence).dg / 1000.0

    @property
    def three_prime_base(self) -> str:
        return self.sequence[-1].upper()


# Secondary-structure thresholds are on free energy, not melting temperature.
# primer3 reports a hairpin Tm that can look alarming for a 3 bp stem whose
# dG is only about -2 kcal/mol, i.e. structures that never actually form at
# annealing temperature. dG is the quantity that decides whether the primer
# folds, so that is what gates the warnings here.
HAIRPIN_DG_LIMIT = -3.0      # kcal/mol
DIMER_DG_LIMIT = -6.0
# How close a dimer's Tm must get to the annealing temperature
# before it is worth warning about.
DIMER_TM_MARGIN = 15.0        # degrees C below the anneal still worth flagging


def check_primer(p: Primer, tm_lo: float = 58.0, tm_hi: float = 68.0) -> Primer:
    """Attach warnings for the usual failure modes. Empty list means clean."""
    w: list[str] = []

    if p.reference:
        # Only structural problems matter for a published sequencing primer.
        if p.hairpin_dg < HAIRPIN_DG_LIMIT:
            w.append(f"hairpin dG {p.hairpin_dg:.1f} kcal/mol is strong")
        p.warnings = w
        return p

    if not (tm_lo <= p.tm_anneal <= tm_hi):
        w.append(f"annealing Tm {p.tm_anneal:.1f} C is outside {tm_lo}-{tm_hi} C")
    if not (40.0 <= p.gc <= 60.0):
        w.append(f"annealing GC {p.gc:.0f}% is outside 40-60%")
    if p.three_prime_base not in "GC":
        w.append("3' end is not G or C (weaker clamp)")
    if p.hairpin_dg < HAIRPIN_DG_LIMIT:
        w.append(f"hairpin dG {p.hairpin_dg:.1f} kcal/mol is strong")
    # A self-dimer only matters if it survives the annealing step. Both
    # restriction sites used for cloning are palindromes (CATATG, CTCGAG),
    # so any primer carrying one can pair with a copy of itself through the
    # site: that is inherent to restriction cloning and no clamp removes it.
    # Judged on dG alone every such primer looks bad. Its melting temperature
    # tells the rest of the story - these dimers melt around 5 C, so at a
    # 60 C anneal they do not exist.
    #
    # This is the same mistake as gating hairpins on Tm, in the opposite
    # direction: one number is not enough. A dimer is flagged when it is both
    # energetically strong and stable near the annealing temperature.
    if p.homodimer_dg < DIMER_DG_LIMIT:
        if p.homodimer_tm >= tm_lo - DIMER_TM_MARGIN:
            w.append(f"self-dimer dG {p.homodimer_dg:.1f} kcal/mol, "
                     f"Tm {p.homodimer_tm:.0f} C, stable near the anneal")
        else:
            p.notes.append(
                f"self-dimer dG {p.homodimer_dg:.1f} kcal/mol but Tm "
                f"{p.homodimer_tm:.0f} C, far below the anneal: this is the "
                "palindromic restriction site pairing with itself and is not "
                "a practical problem")

    seq = p.sequence.upper()
    for base in "ACGT":
        if base * 5 in seq:
            w.append(f"run of 5+ {base}")
            break

    p.warnings = w
    return p


def pair_heterodimer_tm(a: Primer, b: Primer) -> float:
    return primer3.calc_heterodimer(a.sequence, b.sequence).tm


def pair_heterodimer_dg(a: Primer, b: Primer) -> float:
    return primer3.calc_heterodimer(a.sequence, b.sequence).dg / 1000.0


# Clamp candidates: 6 bases 5' of the restriction site so the enzyme has
# something to grip. Several are tried and scored rather than one being
# assumed, because the best choice depends on what it sits next to.
CLAMP_CANDIDATES = (
    "GGAATT", "TATATA", "AAGGAA", "GGTGGT", "TTGGTT", "CGCGGC",
    "AAATAA", "TTAATT", "ACTATA", "AGTACT", "CATCAT", "TGATCA",
)


def _score_primer(p: Primer, target_tm: float) -> float:
    """Lower is better. Mirrors the checks in check_primer, as a cost."""
    cost = 0.0
    cost += 4.0 * abs(p.tm_anneal - target_tm)

    gc = p.gc
    if gc < 40.0:
        cost += 3.0 * (40.0 - gc)
    elif gc > 60.0:
        cost += 3.0 * (gc - 60.0)

    hp = p.hairpin_dg
    if hp < HAIRPIN_DG_LIMIT:
        cost += 15.0 * (HAIRPIN_DG_LIMIT - hp)

    hd = p.homodimer_dg
    if hd < DIMER_DG_LIMIT:
        cost += 25.0 * (DIMER_DG_LIMIT - hd)

    if p.three_prime_base not in "GC":
        cost += 8.0

    cost += 0.5 * max(0, p.length - 30)
    return cost


def _pick_anneal_length(seq: str, target_tm: float, lo: int, hi: int) -> str:
    """Grow the annealing region until Tm is closest to target."""
    best, best_diff = seq[:lo], float("inf")
    for n in range(lo, min(hi, len(seq)) + 1):
        candidate = seq[:n]
        diff = abs(primer3.calc_tm(candidate) - target_tm)
        if diff < best_diff:
            best, best_diff = candidate, diff
    return best


def _search_primer(
    template_5to3: str,
    site: str,
    name: str,
    purpose: str,
    tail_note_fmt: str,
    target_tm: float,
    anneal_lo: int,
    anneal_hi: int,
) -> Primer:
    """
    Grid-search clamp sequence x annealing length, score each candidate,
    return the cheapest. The search space is small enough to enumerate.
    """
    best: Primer | None = None
    best_cost = float("inf")

    for clamp in CLAMP_CANDIDATES:
        for n in range(anneal_lo, min(anneal_hi, len(template_5to3)) + 1):
            anneal = template_5to3[:n]
            cand = Primer(
                name=name,
                sequence=clamp + site + anneal,
                anneal_part=anneal,
                purpose=purpose,
                tail_note=tail_note_fmt.format(clamp=clamp),
            )
            cost = _score_primer(cand, target_tm)
            if cost < best_cost:
                best, best_cost = cand, cost

    assert best is not None
    return check_primer(best)


def design_cloning_primers(
    cds: str,
    label: str = "RsTAL",
    five_site: str = "CATATG",
    five_enzyme: str = "NdeI",
    three_site: str = "CTCGAG",
    three_enzyme: str = "XhoI",
    clamp: str = "GGTGGT",
    target_tm: float = 62.0,
    anneal_lo: int = 18,
    anneal_hi: int = 36,
) -> tuple[Primer, Primer]:
    """
    Build the NdeI/XhoI cloning pair for a CDS that starts with ATG and
    carries NO stop codon (C-terminal fusion).

    Forward: clamp + CATATG + CDS from codon 2 onward.
        The ATG inside the NdeI site IS the start codon, so the primer must
        not repeat it.
    Reverse: clamp + CTCGAG + revcomp(3' end of CDS).
        The insert's last codon runs straight into CTCGAG, which reads
        Leu-Glu and keeps the downstream vector His-tag in frame.
    """
    cds = cds.upper()
    if not cds.startswith("ATG"):
        raise ValueError("CDS must begin with ATG")
    if len(cds) % 3:
        raise ValueError("CDS length is not a multiple of 3")

    fwd = _search_primer(
        template_5to3=cds[3:],
        site=five_site,
        name=f"{label}_NdeI_F",
        purpose=f"5' cloning primer, adds {five_enzyme} site",
        tail_note_fmt=(
            "{clamp} clamp + " + f"{five_site} ({five_enzyme}); "
            "the ATG in the site is the start codon, so codon 1 is not repeated"
        ),
        target_tm=target_tm,
        anneal_lo=anneal_lo,
        anneal_hi=anneal_hi,
    )

    rev = _search_primer(
        template_5to3=revcomp(cds),
        site=three_site,
        name=f"{label}_XhoI_R",
        purpose=f"3' cloning primer, adds {three_enzyme} site, stop codon omitted",
        tail_note_fmt=(
            "{clamp} clamp + " + f"{three_site} ({three_enzyme}); "
            "stop codon deliberately omitted so the vector C-terminal His6 is translated"
        ),
        target_tm=target_tm,
        anneal_lo=anneal_lo,
        anneal_hi=anneal_hi,
    )

    return fwd, rev


def design_internal_seq_primers(
    cds: str,
    label: str = "RsTAL",
    spacing: int = 700,
    target_tm: float = 60.0,
    anneal_lo: int = 18,
    anneal_hi: int = 26,
    search_window: int = 120,
) -> list[Primer]:
    """
    Sanger reads run roughly 700-900 clean bases. The T7 promoter and
    terminator primers cover the two ends, so internal primers are needed
    every ~700 nt to tile the middle of a long insert.

    Rather than taking whatever sits at the exact target position, each
    primer is chosen by searching a window around it. A fixed-position pick
    lands on a GC-rich or self-complementary stretch often enough that it is
    not worth doing.
    """
    cds = cds.upper()
    out: list[Primer] = []
    pos = spacing
    idx = 1
    while pos < len(cds) - anneal_hi:
        best: Primer | None = None
        best_cost = float("inf")
        hi = min(pos + search_window, len(cds) - anneal_lo)
        for start in range(pos, hi):
            for n in range(anneal_lo, anneal_hi + 1):
                if start + n > len(cds):
                    break
                anneal = cds[start:start + n]
                cand = Primer(
                    name=f"{label}_seq_F{idx}",
                    sequence=anneal,
                    anneal_part=anneal,
                    purpose=f"internal sequencing primer, anneals at CDS position {start + 1}",
                )
                cost = _score_primer(cand, target_tm)
                if cost < best_cost:
                    best, best_cost = cand, cost
        if best is not None:
            out.append(check_primer(best))
        pos += spacing
        idx += 1
    return out


def primer_table_rows(primers: list[Primer]) -> list[dict]:
    rows = []
    for p in primers:
        rows.append({
            "Name": p.name,
            "Sequence (5'->3')": p.sequence,
            "Length": p.length,
            "Anneal Tm (C)": round(p.tm_anneal, 1),
            "Full Tm (C)": round(p.tm_full, 1),
            "Anneal GC (%)": round(p.gc, 0),
            "Hairpin dG (kcal/mol)": round(p.hairpin_dg, 1),
            "Self-dimer dG (kcal/mol)": round(p.homodimer_dg, 1),
            "Purpose": p.purpose,
            "Warnings": "; ".join(p.warnings) if p.warnings else "none",
        })
    return rows
