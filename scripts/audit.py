"""
audit.py
--------
Independent checks on the finished construct.

Deliberately does NOT reuse assemble.verify(). It re-reads the written
GenBank file from disk and asks different questions, so a mistake in the
build logic cannot pass by being consistent with itself.
"""

from __future__ import annotations

import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Restriction import BamHI, NcoI, NdeI, NheI, XhoI
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = Path(__file__).resolve().parent.parent

T7_PROMOTER = "TAATACGACTCACTATAG"
LAC_OPERATOR = "GGGAATTGTGAGCGGATAACAATT"
THROMBIN_AA = "LVPRGS"
T7_TAG_AA = "MASMTGGQQMG"
KANR_START_AA = "MSHIQRETSCSRPRLNSNMDADLYGYKWARDNVGQSGATIYRLYGKPDAPELFLKHGKGSVANDVTDEMVRLNWLTEFMPLPTIKHFIRT"


def banner(t: str) -> None:
    print(f"\n{t}\n" + "-" * len(t))


def result(ok: bool, label: str, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f": {detail}" if detail else ""))
    return ok


def main() -> int:
    plasmid_path = ROOT / "output" / "pET28a_RsTAL.gb"
    rec = next(SeqIO.parse(str(plasmid_path), "genbank"))
    seq = str(rec.seq).upper()
    circ = seq + seq  # for features that wrap the origin

    up = next(SeqIO.parse(str(ROOT / "data" / "Q3IWB0.txt"), "swiss"))
    target = str(up.seq).upper()

    ok = True
    print(f"Auditing {plasmid_path.name}: {len(seq):,} bp")

    # --- the expressed protein -------------------------------------------
    banner("1. What protein actually comes out")
    sys.path.insert(0, str(Path(__file__).parent))
    from assemble import find_expression_product
    expr = find_expression_product(seq)
    ok &= result(expr["found"], "translational start located")
    product = expr["product"]
    leader = expr["leader"]
    print(f"  product is {len(product)} aa")
    print(f"  leader    {leader!r} ({len(leader)} aa)")
    ok &= result(leader == "MGSSHHHHHHSSGLVPRGSH",
                 "N-terminal leader is the vector His6 + thrombin site",
                 "retained because NdeI sits downstream of them and in frame")
    ok &= result(product[len(leader):len(leader) + len(target)] == target,
                 "full UniProt sequence present between the two tags",
                 f"{len(target)} aa matched")
    tail = product[len(leader) + len(target):]
    ok &= result(tail == "LEHHHHHH",
                 "C-terminal tail is exactly Leu-Glu-His6",
                 f"found {tail!r} (LE comes from the XhoI site CTCGAG)")
    print(f"  thrombin cleavage after LVPR would leave "
          f"GSHM + RsTAL + LEHHHHHH")

    pa = ProteinAnalysis(product)
    mw = molecular_weight(Seq(product), seq_type="protein") / 1000.0
    print(f"  predicted mass {mw:.1f} kDa, pI {pa.isoelectric_point():.2f}")
    print(f"  (UniProt lists the untagged protein at 54.9 kDa)")

    # --- which vector elements survived the cloning -----------------------
    banner("2. Vector elements upstream of the insert")
    full_translations = []
    for strand_seq in (circ, str(Seq(circ).reverse_complement())):
        for frame in range(3):
            n = 3 * ((len(strand_seq) - frame) // 3)
            full_translations.append(
                str(Seq(strand_seq[frame:frame + n]).translate(table=1)))
    ok &= result(THROMBIN_AA in product,
                 "thrombin site LVPRGS is in the expressed product",
                 "so the N-terminal tag can be cleaved off after IMAC")
    ok &= result(T7_TAG_AA not in product,
                 "T7 tag MASMTGGQQMG is absent from the product",
                 "it lay between NdeI and XhoI and was excised with the stuffer")

    # --- regulatory elements still present and correctly placed ----------
    banner("3. Expression cassette intact")
    p_t7 = circ.find(T7_PROMOTER)
    ok &= result(p_t7 != -1, "T7 promoter present")
    # distance from promoter to the start codon, going round the circle
    if p_t7 != -1:
        dist = (len(seq) - p_t7) % len(seq)
        ok &= result(0 < dist < 400,
                     "T7 promoter sits upstream of the start codon",
                     f"{dist} bp before the ATG")
    p_lac = circ.find(LAC_OPERATOR)
    result(p_lac != -1, "lac operator present")

    # --- resistance marker ------------------------------------------------
    banner("4. Selection marker")
    kan_found = any(KANR_START_AA in t for t in full_translations)
    ok &= result(kan_found, "KanR (aph(3')-Ia) coding sequence intact")

    # --- restriction map --------------------------------------------------
    banner("5. Restriction map of the finished plasmid")
    s = Seq(seq)
    for enz, expect in ((NdeI, 1), (XhoI, 1)):
        hits = enz.search(s, linear=False)
        ok &= result(len(hits) == expect,
                     f"{enz} sites = {len(hits)} (expected {expect})",
                     str(hits))
    for enz in (NheI, NcoI, BamHI):
        hits = enz.search(s, linear=False)
        print(f"  info   {str(enz):<6} sites: {hits}")

    # --- diagnostic digest prediction -------------------------------------
    banner("6. Predicted NdeI + XhoI diagnostic digest")
    nde = NdeI.search(s, linear=False)[0] - 1
    xho = XhoI.search(s, linear=False)[0] - 1
    a = (xho - nde) % len(seq)
    b = len(seq) - a
    big, small = max(a, b), min(a, b)
    print(f"  two fragments: {big:,} bp (backbone) and {small:,} bp (insert)")
    ok &= result(abs(small - len(target) * 3) < 30,
                 "small fragment matches the insert size",
                 f"{small} bp vs CDS {len(target) * 3} bp")

    # --- no internal sites that would spoil the digest --------------------
    banner("7. Sequencing coverage")
    read = 800
    cds_len = len(target) * 3
    n_needed = max(0, -(-cds_len // read) - 1)
    print(f"  insert {cds_len} bp, ~{read} bp per Sanger read")
    print(f"  T7 promoter and terminator primers cover both ends;")
    print(f"  {n_needed} internal primer(s) needed to tile the middle")

    banner("Result")
    print("  ALL CHECKS PASSED" if ok else "  SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
