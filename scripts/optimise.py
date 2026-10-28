"""
optimise.py
-----------
Reverse-translate and codon-optimise a protein for E. coli expression,
while simultaneously forbidding the restriction sites needed for cloning.

Doing both in one solver matters. If you optimise first and strip sites
afterwards you can reintroduce a rare codon; if you strip sites first and
then optimise you can recreate the site. DNA Chisel treats the forbidden
sites as hard constraints and codon usage as an objective, so the result
satisfies both at once.

Hard constraints
    - translation must be identical to the input protein
    - no occurrence of any forbidden restriction site (both strands)
    - GC between 30% and 70% in every 50 nt window
    - no stable hairpins

Objective
    - match E. coli codon usage
"""

from __future__ import annotations

import numpy as np

from dnachisel import (
    AvoidHairpins,
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceTranslation,
    reverse_translate,
)

# Sites we must not have inside the insert.
# NdeI and XhoI are the cloning sites; the rest are the other pET-28a MCS
# enzymes, kept clear so the construct stays easy to work with later.
FORBIDDEN_SITES = {
    "NdeI": "CATATG",
    "XhoI": "CTCGAG",
    "NcoI": "CCATGG",
    "NheI": "GCTAGC",
    "BamHI": "GGATCC",
    "EcoRI": "GAATTC",
    "SacI": "GAGCTC",
    "SalI": "GTCGAC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "EcoRV": "GATATC",
}


def optimise_cds(
    protein: str,
    forbidden: dict[str, str] | None = None,
    gc_min: float = 0.30,
    gc_max: float = 0.70,
    gc_window: int = 50,
    global_gc: tuple[float, float] = (0.35, 0.62),
    ramp_codons: int = 15,
    ramp_gc: tuple[float, float] = (0.25, 0.50),
    seed: int = 1,
) -> dict:
    """
    Returns a dict with the naive reverse translation, the optimised CDS,
    and the solver's constraint/objective reports.

    The returned CDS has no stop codon: the construct is a C-terminal fusion,
    so the stop is supplied by the vector downstream of the insert.

    `ramp_codons` / `ramp_gc` implement the 5'-end ramp. The first ~15 codons
    of a bacterial CDS should be AT-biased and structure-poor: a stable
    hairpin over the ribosome binding site and start codon blocks initiation,
    and that single effect sinks more expression constructs than poor codon
    usage across the rest of the gene ever does. Naively maximising CAI makes
    this worse for a protein from a high-GC organism, because the preferred
    codons are themselves GC-rich.
    """
    if forbidden is None:
        forbidden = FORBIDDEN_SITES

    protein = protein.strip().upper()
    if "*" in protein:
        raise ValueError("Strip the stop (*) from the protein before optimising.")
    if not protein.startswith("M"):
        raise ValueError("Protein does not start with Met.")

    # DNA Chisel's solver is stochastic. Without a fixed seed the same input
    # gives a different (still valid) sequence every run, which would make the
    # reported CAI, the primers and the final plasmid irreproducible.
    np.random.seed(seed)

    naive = reverse_translate(protein)

    constraints = [EnforceTranslation()]
    for name, site in forbidden.items():
        constraints.append(AvoidPattern(site))
    constraints.append(
        EnforceGCContent(mini=gc_min, maxi=gc_max, window=gc_window)
    )
    constraints.append(
        EnforceGCContent(mini=global_gc[0], maxi=global_gc[1])
    )
    # 5' translation-initiation ramp
    constraints.append(
        EnforceGCContent(
            mini=ramp_gc[0], maxi=ramp_gc[1],
            location=(0, ramp_codons * 3),
        )
    )
    constraints.append(AvoidHairpins(stem_size=10, hairpin_window=60))

    # The two ends of the CDS are also the primer-binding sites, so they get
    # extra treatment: no short hairpins (which would fold the primer itself)
    # and a GC content a primer can actually be designed from.
    n = len(naive)
    constraints.append(AvoidHairpins(stem_size=6, hairpin_window=40, location=(0, 60)))
    constraints.append(AvoidHairpins(stem_size=6, hairpin_window=40, location=(n - 60, n)))
    constraints.append(EnforceGCContent(mini=0.45, maxi=0.62, location=(n - 45, n)))

    problem = DnaOptimizationProblem(
        sequence=naive,
        constraints=constraints,
        objectives=[CodonOptimize(species="e_coli")],
        logger=None,
    )

    problem.resolve_constraints()
    problem.optimize()

    return {
        "protein": protein,
        "naive_cds": naive,
        "optimised_cds": problem.sequence,
        "constraints_report": problem.constraints_text_summary(),
        "objectives_report": problem.objectives_text_summary(),
        "all_constraints_pass": problem.all_constraints_pass(),
    }


def find_sites(seq: str, sites: dict[str, str] | None = None) -> dict[str, list[int]]:
    """
    Scan both strands for each site. Returns {enzyme: [1-based positions]}.
    Only enzymes that are actually present appear in the result.
    """
    if sites is None:
        sites = FORBIDDEN_SITES

    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    seq = seq.upper()
    rc = seq.translate(complement)[::-1]

    hits: dict[str, list[int]] = {}
    for name, site in sites.items():
        positions: list[int] = []
        start = seq.find(site)
        while start != -1:
            positions.append(start + 1)
            start = seq.find(site, start + 1)
        # Palindromic sites would be double-counted, so only scan the reverse
        # strand for non-palindromic recognition sequences.
        if site != site.translate(complement)[::-1]:
            start = rc.find(site)
            while start != -1:
                positions.append(-(len(seq) - start - len(site) + 1))
                start = rc.find(site, start + 1)
        if positions:
            hits[name] = sorted(positions, key=abs)
    return hits
