"""Script tests for pre_irma.py, which is used to process input files before the main IRMA pipeline run."""
from pathlib import Path

import pytest

from pre_irma import (
    get_species_and_seq_type,
    group_by_sample,
    has_r1_r2,
    read_number,
    sample_name,
)


@pytest.mark.parametrize("name, expected", [
    ("26-040001_S1_L001_R1_001.fastq.gz", 1),   # MiSeq / bcl2fastq
    ("26-040001_S1_L001_R2_001.fastq.gz", 2),
    ("sample_R1.fastq.gz", 1),
    ("sample_R2.fq.gz", 2),
    ("sample_1.fastq", 1),
    ("sample_2.fq", 2),
    ("sample.fastq.gz", None),
    ("26-040002_S2_L001_I1_001.fastq.gz", None),  # index read
])
def test_read_number(name, expected):
    assert read_number(name) == expected


def test_lone_r1_is_not_a_pair():
    # Old matcher: 'r2' substring / "-2." anywhere counted -> this looked paired.
    assert has_r1_r2(["ctr2-2_S1_L001_R1_001.fastq.gz"]) == (True, False)


def test_miseq_pair():
    files = ["x/26-040001_S1_L001_R1_001.fastq.gz", "x/26-040001_S1_L001_R2_001.fastq.gz"]
    assert has_r1_r2(files) == (True, True)


def test_miseq_sample_names_and_grouping():
    files = [Path(f"/in/influenza/miseq/run/Data/Intensities/BaseCalls/{s}_L001_R{r}_001.fastq.gz")
             for s in ("26-040001_S1", "26-040002_S2", "Undetermined_S0") for r in (1, 2)]
    groups = group_by_sample(files)
    assert sorted(groups) == ["26-040001", "26-040002", "Undetermined"]
    assert all(len(v) == 2 for v in groups.values())
    assert sample_name("26-040001_S1_L001_R1_001.fastq.gz") == "26-040001"


def test_species_and_seq_type_for_miseq_run_folder(tmp_path: Path):
    f = tmp_path / "influenza" / "MiSeq" / "run" / "Data" / "x_S1_L001_R1_001.fastq.gz"
    f.parent.mkdir(parents=True)
    f.write_text("")
    assert get_species_and_seq_type([str(f)], tmp_path) == ("FLU", "miseq")


def test_combine_fastqs_mixes_gzip_and_plain(tmp_path: Path):
    import gzip

    from pre_irma import combine_fastqs
    with gzip.open(tmp_path / "a.fastq.gz", "wt") as fh:
        fh.write("@a\nAC\n+\nII\n")
    (tmp_path / "b.fastq").write_text("@b\nGT\n+\nII\n")
    out = combine_fastqs([tmp_path / "a.fastq.gz", tmp_path / "b.fastq"], tmp_path / "out" / "x.fastq.gz")
    with gzip.open(out, "rt") as fh:
        assert fh.read() == "@a\nAC\n+\nII\n@b\nGT\n+\nII\n"
    assert not (tmp_path / "out" / "x.fastq.gz.part").exists()


def test_unstable_files_flags_empty_and_missing(tmp_path: Path):
    from pre_irma import unstable_files
    good, empty = tmp_path / "good.fq", tmp_path / "empty.fq"
    good.write_text("@r\nA\n+\nI\n")
    empty.write_text("")
    assert unstable_files([good, empty, tmp_path / "gone.fq"], wait_seconds=0) == [empty, tmp_path / "gone.fq"]


def test_minion_barcodes_and_unclassified_grouping():
    run = Path("/in/influenza/minion/RUN")
    files = [run / "fastq_pass/barcode01/a_0.fastq.gz", run / "fastq_fail/barcode01/a_0.fastq.gz",
             run / "fastq_pass/unclassified/u_0.fastq.gz", run / "fastq_pass/unclassified/u_1.fastq.gz"]
    groups = group_by_sample(files)
    assert sorted(groups) == ["barcode01", "unclassified"]
    assert len(groups["barcode01"]) == 2
