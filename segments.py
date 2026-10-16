"""Recognise influenza segments in IRMA consensus file names.

Every caller goes through ``consensus_segments()`` (or ``consensus_fastas()``
and ``segment_from_filename()``), so the naming rules live in one place.
``read_fasta()`` reads the consensus files themselves.

How IRMA names ``amended_consensus/`` files (``IRMA_RES/scripts/
amendConsensus.pl``, https://github.com/CDCgov/irma): each ``SEG_NUMBERS``
``key:value`` pair is tried in order, the key is matched as a regex against
the gene name (``A_HA_H3``, ``B_PB1``, ...), and the first match appends
``_<value>`` to the sample name. Modules without ``SEG_NUMBERS`` (CoV, RSV)
write a plain ``<sample>.fa``.

* The lab's ``SEG_NUMBERS`` (``B_PB1:B_PB1,...,HA:HA,NA:NA,MP:MP,NS:NS``)
  gives names: ``26-040001_HA.fa``, ``26-040001_A_PB1.fa``. Tokens are
  matched exactly, so the ``A_``/``B_`` prefix is fine.
* IRMA's stock FLU setting (``IRMA_RES/modules/FLU/init.sh``) gives numbers:
  ``B_PB1:1,B_PB2:2,A_PB2:1,A_PB1:2,PA:3,HA:4,NP:5,NA:6,M:7,NS:8``.
  Numbers 1 and 2 mean different segments for influenza A and B, so they
  can only be decoded when the type is known (see ``flu_type_of()``).
* With ``PADDED_CONSENSUS=1``, every ``<x>.fa`` also gets an N-padded
  ``<x>.pad.fa``. Counting those as consensus would double the segments
  and inflate N%.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

logger = logging.getLogger(__name__)

FLU_SEGMENTS = ("HA", "NA", "PB1", "PB2", "PA", "MP", "NS", "NP")

# IRMA's stock SEG_NUMBERS for the FLU module. 1 and 2 depend on the type.
IRMA_SEGMENT_NUMBERS = {3: "PA", 4: "HA", 5: "NP", 6: "NA", 7: "MP", 8: "NS"}
IRMA_POLYMERASE_NUMBERS = {
    "A": {1: "PB2", 2: "PB1"},
    "B": {1: "PB1", 2: "PB2"},
}

FASTA_SUFFIXES = (".fa", ".fasta", ".fas")
PADDED_SUFFIX = ".pad.fa"


def is_padded(path: Path | str) -> bool:
    """Check whether a file is one of IRMA's N-padded consensus copies.

    Args:
        path: File name or path.

    Returns:
        ``True`` for ``*.pad.fa`` files.
    """
    return Path(path).name.lower().endswith(PADDED_SUFFIX)


def consensus_fastas(amended_dir: Path) -> list[Path]:
    """List the real (non-padded) consensus FASTA files in a folder.

    Args:
        amended_dir: IRMA's ``amended_consensus`` folder.

    Returns:
        The FASTA files, sorted; an empty list if the folder does not exist.
    """
    if not amended_dir.is_dir():
        return []
    return sorted(
        p for p in amended_dir.iterdir()
        if p.is_file() and p.name.lower().endswith(FASTA_SUFFIXES) and not is_padded(p)
    )



def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read the records of a FASTA file.

    Args:
        path: FASTA file.

    Returns:
        ``(header, sequence)`` per record, in file order. The header is
        without the leading ``>``; the sequence lines are joined.
    """
    records: list[tuple[str, list[str]]] = []
    with open(path) as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                records.append((line[1:], []))
            elif line and records:
                records[-1][1].append(line)
    return [(header, "".join(lines)) for header, lines in records]

def _strip_fasta_suffix(name: str) -> str:
    """Remove a ``.pad.fa``, ``.fa``, ``.fasta`` or ``.fas`` suffix.

    Args:
        name: File name.

    Returns:
        ``name`` without the suffix (matched case-insensitively), or
        unchanged if it has none.
    """
    lower = name.lower()
    for suffix in (PADDED_SUFFIX, *FASTA_SUFFIXES):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def segment_from_filename(name: str, flu_type: str | None = None) -> str | None:
    """Return the segment (HA, NA, PB1, ...) a consensus file name refers to.

    Handles both IRMA naming styles::

        A_HA_H3.fasta, sample_HA.fa, sample_B_PB1.fa   -> named segment
        26-040001_4.fa, 26-040001_4.pad.fa             -> IRMA segment number

    The name is split on ``_``, ``-`` and ``.`` and tokens are matched
    exactly, so a sample name that merely contains "NA" or "HA" is not
    mistaken for a segment.

    Args:
        name: File name or path.
        flu_type: ``"A"`` or ``"B"``. Needed to decode numbers 1 and 2
            (polymerase genes); without it those return ``None``.

    Returns:
        Segment name, or ``None`` if it cannot be determined.
    """
    tokens = [t.upper() for t in re.split(r"[_.\-]", _strip_fasta_suffix(Path(name).name)) if t]
    for token in tokens:
        if token in FLU_SEGMENTS:
            return token
    if tokens and tokens[-1].isdigit():
        number = int(tokens[-1])
        if number in IRMA_SEGMENT_NUMBERS:
            return IRMA_SEGMENT_NUMBERS[number]
        return IRMA_POLYMERASE_NUMBERS.get((flu_type or "").upper(), {}).get(number)
    return None


def flu_type_of(irma_sample_dir: Path) -> str | None:
    """Find the influenza type of an IRMA sample folder.

    The type is read from IRMA's per-gene files in the folder
    (``A_HA_H3.fasta``, ``B_NA.fasta``, ...).

    Args:
        irma_sample_dir: IRMA output folder of one sample.

    Returns:
        ``"A"`` or ``"B"``, or ``None`` if there are no per-gene files.
    """
    for fasta in sorted(irma_sample_dir.glob("*.fasta")):
        first = fasta.name.split("_", 1)[0].upper()
        if first in ("A", "B"):
            return first
    return None


def consensus_segments(amended_dir: Path) -> list[tuple[Path, str]]:
    """Pair each real consensus file in a folder with its segment.

    The influenza type is read from the parent IRMA sample folder, so
    numbered names can be decoded. Files whose segment cannot be determined
    are logged and skipped.

    Args:
        amended_dir: IRMA's ``amended_consensus`` folder.

    Returns:
        ``(fasta, segment)`` for each file whose segment is known.
    """
    flu_type = flu_type_of(amended_dir.parent)
    result = []
    for fasta in consensus_fastas(amended_dir):
        seg = segment_from_filename(fasta.name, flu_type)
        if seg is None:
            logger.warning("Cannot tell which influenza segment %s is (type %s); skipping it.",
                           fasta, flu_type or "unknown")
            continue
        result.append((fasta, seg))
    return result