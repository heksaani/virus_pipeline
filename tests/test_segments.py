"""CHANGE (2026-10-06): NEW. Tests for segments.py (IRMA filename parsing)."""
from pathlib import Path

import pytest

from segments import (
    consensus_fastas,
    consensus_segments,
    flu_type_of,
    is_padded,
    segment_from_filename,
)


@pytest.mark.parametrize("name, expected", [
    # IRMA amended_consensus: numbered by SEG_NUMBERS
    ("26-040001_3.fa", "PA"),
    ("26-040001_4.fa", "HA"),
    ("26-040001_5.fa", "NP"),
    ("26-040001_6.fa", "NA"),
    ("26-040001_7.fa", "MP"),
    ("26-040001_8.fa", "NS"),
    ("26-040001_4.pad.fa", "HA"),
    # lab SEG_NUMBERS (flu-amd/IRMA_RES/modules/FLU/init.sh, 2026-10)
    ("26-040001_HA.fa", "HA"),
    ("26-040001_NA.pad.fa", "NA"),
    ("26-040001_A_PB1.fa", "PB1"),
    ("26-040001_A_PB2.fa", "PB2"),
    ("26-040001_B_PB1.fa", "PB1"),
    ("26-040001_B_PB2.pad.fa", "PB2"),
    # IRMA sample root: named
    ("A_HA_H3.fasta", "HA"),
    ("A_NA_N2.fasta", "NA"),
    ("B_PB1.fasta", "PB1"),
    ("sample_MP.fa", "MP"),
])
def test_segment_from_filename(name, expected):
    assert segment_from_filename(name) == expected


@pytest.mark.parametrize("name", [
    "26-040001.fa",        # COVID single genome
    "26-040001_9.fa",      # not a flu segment number
    "HANNA_report.fa",     # contains "HA"/"NA" but not as a token
    "notes.txt",
])
def test_segment_from_filename_unknown(name):
    assert segment_from_filename(name) is None


def test_consensus_fastas_skips_padded_and_non_fasta(tmp_path: Path):
    for name in ["s_4.fa", "s_4.pad.fa", "s_6.fa", "s_6.pad.fa", "s_4.a2m", "README"]:
        (tmp_path / name).write_text(">x\nACGT\n")
    assert [p.name for p in consensus_fastas(tmp_path)] == ["s_4.fa", "s_6.fa"]


def test_consensus_fastas_missing_dir(tmp_path: Path):
    assert consensus_fastas(tmp_path / "nope") == []


def test_is_padded():
    assert is_padded("x_4.pad.fa")
    assert not is_padded("x_4.fa")


# IRMA's stock SEG_NUMBERS (github.com/CDCgov/irma, IRMA_RES/modules/FLU/init.sh):
# B_PB1:1,B_PB2:2,A_PB2:1,A_PB1:2 -- numbers 1 and 2 depend on the type.
@pytest.mark.parametrize("name, flu_type, expected", [
    ("s_1.fa", "A", "PB2"), ("s_2.fa", "A", "PB1"),
    ("s_1.fa", "B", "PB1"), ("s_2.fa", "B", "PB2"),
    ("s_1.fa", None, None),       # never guess
    ("s_4.fa", None, "HA"),       # 3-8 are the same for A and B
])
def test_numbered_polymerase_segments_need_type(name, flu_type, expected):
    assert segment_from_filename(name, flu_type) == expected


def test_consensus_segments_reads_type_from_irma_folder(tmp_path: Path):
    irma = tmp_path / "s"
    (irma / "amended_consensus").mkdir(parents=True)
    (irma / "B_HA.fasta").write_text(">B_HA\nACGT\n")
    for n in (1, 2, 4):
        (irma / "amended_consensus" / f"s_{n}.fa").write_text(">x\nACGT\n")
    assert flu_type_of(irma) == "B"
    assert [seg for _, seg in consensus_segments(irma / "amended_consensus")] == ["PB1", "PB2", "HA"]
