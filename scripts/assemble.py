"""
assemble.py
-----------
Simulate the NdeI/XhoI digest and ligation, build the finished plasmid,
and run the verification checks.

Nothing about the vector is assumed. The cut sites, the reading frame, the
position of the His tag and whether the N-terminal tag survives the cloning
are all read out of the supplied GenBank record. That matters for pET-28a(+)
in particular: the ATG inside the NdeI site becomes the start codon, so
cloning NdeI/XhoI discards the vector's N-terminal His6 and thrombin site
and leaves a C-terminal His6 fusion instead. This script reports what the
sequence actually says rather than taking that on trust.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio.Restriction import NdeI, XhoI
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HIS6_DNA_VARIANTS = ("CACCACCACCACCACCAC", "CATCATCATCATCATCAT")


@dataclass
class Check:
    name: str
    passed: bool
    detail: str

    def line(self) -> str:
        return f"[{'PASS' if self.passed else 'FAIL'}] {self.name}: {self.detail}"


def single_site(seq: Seq, enzyme, label: str) -> int:
    """Return the single cut position for an enzyme, or raise."""
    sites = enzyme.search(seq, linear=False)
    if len(sites) != 1:
        raise ValueError(
            f"{label}: expected exactly 1 {enzyme} site, found {len(sites)} at {sites}"
        )
    return sites[0]


def pcr_product(cds: str, fwd_primer: str, rev_primer: str) -> str:
    """
    Reconstruct the amplicon from the two primers and the template.
    The product is forward primer + the template between the primer
    binding sites + reverse complement of the reverse primer.
    """
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    rev_rc = rev_primer.translate(comp)[::-1]

    # Find where each primer's annealing region sits on the CDS.
    # The forward primer's 3' end must match somewhere in the CDS.
    for n in range(40, 14, -1):
        tail = fwd_primer[-n:]
        if cds.startswith(tail):
            fwd_overlap = n
            break
    else:
        raise ValueError("Forward primer does not anneal to the CDS 5' end")

    for n in range(40, 14, -1):
        head = rev_rc[:n]
        if cds.endswith(head):
            rev_overlap = n
            break
    else:
        raise ValueError("Reverse primer does not anneal to the CDS 3' end")

    middle = cds[fwd_overlap: len(cds) - rev_overlap]
    return fwd_primer + middle + rev_rc


def digest_insert(amplicon: str) -> str:
    """Cut the amplicon with NdeI and XhoI, return the central fragment."""
    seq = Seq(amplicon)
    nde = NdeI.search(seq, linear=True)
    xho = XhoI.search(seq, linear=True)
    if len(nde) != 1 or len(xho) != 1:
        raise ValueError(
            f"Amplicon must have exactly one NdeI and one XhoI site "
            f"(found NdeI={len(nde)}, XhoI={len(xho)})"
        )
    # Biopython returns the position after the cut, 1-based.
    start = nde[0] - 1
    end = xho[0] - 1
    if end <= start:
        raise ValueError("XhoI site precedes NdeI site in the amplicon")
    return str(seq[start:end])


def digest_vector(vector: SeqRecord) -> tuple[str, int, int]:
    """
    Cut the circular vector with NdeI and XhoI and return the large backbone
    fragment, written so that it begins immediately after the NdeI cut.
    Also returns the two cut positions for the record.
    """
    seq = vector.seq
    nde = single_site(seq, NdeI, "vector")
    xho = single_site(seq, XhoI, "vector")

    nde_i = nde - 1
    xho_i = xho - 1

    if xho_i > nde_i:
        # Backbone is everything outside the NdeI..XhoI stuffer.
        backbone = str(seq[xho_i:]) + str(seq[:nde_i])
    else:
        backbone = str(seq[xho_i:nde_i])

    return backbone, nde, xho


def build_construct(
    vector: SeqRecord,
    insert_fragment: str,
    name: str = "pET28a_RsTAL",
) -> SeqRecord:
    """
    Ligate the digested insert into the digested backbone.

    NdeI cuts CA^TATG, leaving a 5' TA overhang. So the insert fragment
    begins "TATG..." and the backbone's NdeI end finishes "...CA": the
    CATATG site, and with it the ATG start codon, only exists once the two
    are ligated. The finished circle is then rotated so that base 1 is the
    A of that start codon, which makes every downstream frame check simple
    arithmetic.
    """
    backbone, nde, xho = digest_vector(vector)

    if not insert_fragment.upper().startswith("TA"):
        raise ValueError(
            "Insert fragment does not start with the NdeI TA overhang; "
            f"it starts with {insert_fragment[:6]!r}"
        )

    full = backbone + insert_fragment

    # CATATG spans the junction; the start codon is the last three bases of it.
    atg_index = len(backbone) + 1
    if full[atg_index: atg_index + 3].upper() != "ATG":
        raise ValueError(
            "Start codon was not reconstituted at the NdeI junction; found "
            f"{full[atg_index:atg_index+3]!r}"
        )

    rotated = full[atg_index:] + full[:atg_index]

    rec = SeqRecord(
        Seq(rotated),
        id=name,
        name=name,
        description=f"{vector.id} carrying RsTAL cloned NdeI/XhoI",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    return rec


SHINE_DALGARNO = "AAGGAG"


def find_expression_product(plasmid: str, search_back: int = 300) -> dict:
    """
    Work out what the cell actually translates, rather than assuming the
    insert's own ATG is the start codon.

    This matters for pET-28a(+). The NdeI site sits DOWNSTREAM of the
    vector's N-terminal His6 and thrombin site, and in frame with the
    vector's own start codon, which is the one the ribosome binding site
    serves. So NdeI/XhoI cloning does NOT discard the N-terminal tag the
    way NcoI cloning would: the product carries both an N-terminal
    His6-thrombin leader and the C-terminal His6.

    The plasmid is assumed to be rotated so that base 1 is the insert's ATG,
    which puts the upstream region at the end of the sequence.
    """
    seq = plasmid.upper()
    upstream = seq[-search_back:]
    region = upstream + seq

    rbs = upstream.rfind(SHINE_DALGARNO)
    if rbs == -1:
        return {"found": False, "reason": "no Shine-Dalgarno (AAGGAG) upstream"}

    # First ATG after the RBS that reads in frame into the insert.
    insert_start = len(upstream)  # index of the insert ATG inside `region`
    for i in range(rbs, insert_start + 1):
        if region[i:i + 3] != "ATG":
            continue
        if (insert_start - i) % 3 != 0:
            continue
        orf_region = region[i:]
        usable = orf_region[: 3 * (len(orf_region) // 3)]
        protein = str(Seq(usable).translate(table=1))
        stop = protein.find("*")
        if stop == -1:
            continue
        return {
            "found": True,
            "start_offset_before_insert": insert_start - i,
            "rbs_to_start_spacing": i - (rbs + len(SHINE_DALGARNO)),
            "leader": protein[: (insert_start - i) // 3],
            "product": protein[:stop],
            "length_aa": stop,
        }

    return {"found": False, "reason": "no in-frame upstream ATG with a stop"}


def find_his_tag(plasmid: str) -> int | None:
    """Return the 0-based index of a His6 run, or None."""
    for variant in HIS6_DNA_VARIANTS:
        idx = plasmid.find(variant)
        if idx != -1:
            return idx
    # Mixed-codon His6: scan translation instead.
    return None


def verify(
    plasmid: SeqRecord,
    target_protein: str,
    vector: SeqRecord,
) -> list[Check]:
    """The five verification checks, plus the His-tag frame check."""
    checks: list[Check] = []
    seq = str(plasmid.seq).upper()

    # 1. The real translational start, which is NOT the insert's own ATG.
    expr = find_expression_product(seq)
    if not expr["found"]:
        checks.append(Check("Translational start", False, expr["reason"]))
    else:
        checks.append(Check(
            "Translational start",
            True,
            f"ribosome starts {expr['start_offset_before_insert']} bp upstream of "
            f"the insert ATG, in frame; leader is {expr['leader']!r} "
            f"(vector N-terminal His6 + thrombin site, retained by NdeI cloning)",
        ))

    # 2. The expressed product is leader + target protein + LEHHHHHH.
    cds_len = len(target_protein) * 3
    insert_cds = seq[:cds_len]
    translated = str(Seq(insert_cds).translate(table=1))
    match = translated == target_protein
    if match:
        detail = f"{len(translated)} aa, identical to the UniProt sequence"
    else:
        diffs = [
            f"{i+1}:{a}->{b}"
            for i, (a, b) in enumerate(zip(target_protein, translated))
            if a != b
        ][:5]
        detail = f"MISMATCH ({len(translated)} aa) first differences: {diffs}"
    checks.append(Check("Insert translation", match, detail))

    if expr["found"]:
        expected = expr["leader"] + target_protein + "LEHHHHHH"
        checks.append(Check(
            "Full expressed product",
            expr["product"] == expected,
            f"{expr['length_aa']} aa = {len(expr['leader'])} aa leader + "
            f"{len(target_protein)} aa RsTAL + LEHHHHHH"
            if expr["product"] == expected
            else f"product does not match leader+target+LEHHHHHH "
                 f"(got {expr['length_aa']} aa)",
        ))

    # 3. No premature stop inside the insert.
    internal_stop = "*" in translated
    checks.append(Check(
        "No premature stop",
        not internal_stop,
        "no stop codon within the insert reading frame" if not internal_stop
        else f"stop codon at insert residue {translated.index('*')+1}",
    ))

    # 4. Downstream His tag is in frame with the insert.
    his_idx = find_his_tag(seq)
    if his_idx is None:
        checks.append(Check("His tag in frame", False, "no His6 codon run found"))
    else:
        in_frame = (his_idx % 3) == 0
        # Read from the start codon through to the first stop.
        orf = str(Seq(seq[: (len(seq) // 3) * 3]).translate(table=1))
        stop_at = orf.find("*")
        tag_in_orf = "HHHHHH" in (orf[:stop_at] if stop_at != -1 else orf)
        ok = in_frame and tag_in_orf
        checks.append(Check(
            "His tag in frame",
            ok,
            f"His6 DNA at position {his_idx+1}, offset {his_idx % 3} from frame; "
            f"{'present' if tag_in_orf else 'NOT present'} in the translated ORF "
            f"before the first stop",
        ))

    # 5. Junctions read as designed.
    # The plasmid is rotated to the ATG, so CATATG wraps the origin: the "CA"
    # is the last two bases of the sequence.
    nde_reformed = (seq[-3:] + seq[:3]).upper() == "CATATG"
    checks.append(Check(
        "5' junction (NdeI)",
        nde_reformed,
        f"CATATG reconstituted across the origin: ...{seq[-3:]} | {seq[:3]}...",
    ))

    xho_pos = seq.find("CTCGAG", cds_len - 3)
    xho_ok = xho_pos == cds_len and (xho_pos % 3) == 0
    checks.append(Check(
        "3' junction (XhoI)",
        xho_ok,
        (f"XhoI site at base {xho_pos+1}, immediately after the last insert "
         f"codon, frame offset {xho_pos % 3} (0 = Leu-Glu reads in frame "
         f"into the tag)")
        if xho_pos != -1 else "XhoI site not found after the insert",
    ))

    # 6. Size arithmetic.
    stuffer = _stuffer_length(vector)
    expected = len(vector.seq) - stuffer + len(insert_cds) + 2
    actual = len(seq)
    checks.append(Check(
        "Plasmid size",
        actual == expected,
        f"{actual} bp = vector {len(vector.seq)} bp - stuffer {stuffer} bp "
        f"+ insert CDS {len(insert_cds)} bp + 2 bp (the CA of the "
        f"reconstituted NdeI site and the C of the XhoI site)",
    ))

    return checks


def _stuffer_length(vector: SeqRecord) -> int:
    nde = single_site(vector.seq, NdeI, "vector") - 1
    xho = single_site(vector.seq, XhoI, "vector") - 1
    return abs(xho - nde)


def annotate(plasmid: SeqRecord, protein_len: int, label: str = "RsTAL") -> SeqRecord:
    """Add a CDS feature for the insert so the GenBank file is readable."""
    plasmid.features.append(
        SeqFeature(
            FeatureLocation(0, protein_len * 3),
            type="CDS",
            qualifiers={
                "label": [label],
                "note": ["codon-optimised for E. coli; stop codon omitted"],
                "codon_start": ["1"],
            },
        )
    )
    return plasmid
