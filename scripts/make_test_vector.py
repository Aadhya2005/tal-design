"""
make_test_vector.py
-------------------
Builds a SYNTHETIC TEST FIXTURE that mimics the parts of pET-28a(+) the
pipeline touches: a single NdeI site, a short stuffer, a single XhoI site,
a C-terminal His6 and a stop codon, inside a circular backbone.

This is NOT pET-28a(+). It exists only so the assembly and verification code
can be proven correct before the real vector arrives. The real run uses the
GenBank file downloaded from SnapGene or Addgene.
"""

from __future__ import annotations

import random

from Bio.Restriction import NdeI, XhoI
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

FORBIDDEN = ("CATATG", "CTCGAG")


def clean_filler(n: int, rng: random.Random) -> str:
    """Random filler guaranteed to contain neither cloning site."""
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(n))
        if not any(f in s for f in FORBIDDEN):
            return s


def build(seed: int = 42) -> SeqRecord:
    rng = random.Random(seed)

    upstream = clean_filler(600, rng)          # stands in for ori / KanR
    t7 = "TAATACGACTCACTATAGGG"                # T7 promoter primer site
    lac_rbs = clean_filler(40, rng) + "AAGGAG" + "ATATA"
    nde = "CATATG"
    stuffer = clean_filler(99, rng)            # multiple of 3 keeps frame tidy
    xho = "CTCGAG"
    his_stop = "CACCACCACCACCACCAC" + "TGA"
    t7_term = "GCTAGTTATTGCTCAGCGG"
    downstream = clean_filler(500, rng)

    seq = upstream + t7 + lac_rbs + nde + stuffer + xho + his_stop + t7_term + downstream

    rec = SeqRecord(
        Seq(seq),
        id="TEST_FIXTURE_vector",
        name="TEST_FIXTURE",
        description="SYNTHETIC TEST FIXTURE - not pET-28a(+)",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )

    assert len(NdeI.search(rec.seq, linear=False)) == 1, "fixture needs one NdeI site"
    assert len(XhoI.search(rec.seq, linear=False)) == 1, "fixture needs one XhoI site"
    return rec


if __name__ == "__main__":
    v = build()
    print(f"fixture length : {len(v.seq)} bp")
    print(f"NdeI cut (1-based): {NdeI.search(v.seq, linear=False)}")
    print(f"XhoI cut (1-based): {XhoI.search(v.seq, linear=False)}")
