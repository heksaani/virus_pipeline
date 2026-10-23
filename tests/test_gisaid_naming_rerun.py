"""Some tests for gisaid_naming_rerun.py, which is used to name samples in a rerun of an IRMA run."""
from pathlib import Path

import pytest

from gisaid_naming import is_thl_sample_name, process_sample_folder
from gisaid_naming_rerun import find_irma_sample_dirs, is_minion_barcode


def test_rerun_all_skips_minion_barcodes(tmp_path: Path):
    for sample in ("26-040001", "barcode01", "Barcode02"):
        (tmp_path / sample / "irma_results" / sample / "amended_consensus").mkdir(parents=True)
    assert [name for name, _ in find_irma_sample_dirs(tmp_path)] == ["26-040001"]


def test_is_minion_barcode():
    assert is_minion_barcode("barcode12")
    assert not is_minion_barcode("26-040001")

@pytest.mark.parametrize("name, expected", [
    ("26-057848", True), ("2026-1234", True),
    ("Neg-KO-IAV", False), ("barcode01", False), ("Undetermined", False), ("26-057848b", False),
    # CHANGE (2026-10-23): VI suffix accepted, nothing else after the digits.
    ("2026-00418VI", True), ("26-057848VI", True),
    ("2026-00418vi", False), ("2026-00418VIX", False), ("2026-VI", False),
])
def test_is_thl_sample_name(name, expected):
    assert is_thl_sample_name(name) is expected


def test_negative_control_is_not_named_and_central_list_untouched(tmp_path: Path):
    irma = tmp_path / "Neg-KO-IAV"
    (irma / "amended_consensus").mkdir(parents=True)
    (irma / "A_HA_H3.fasta").write_text(">A_HA_H3\nACGT\n")      # contamination: HA + NA present
    (irma / "A_NA_N2.fasta").write_text(">A_NA_N2\nACGT\n")
    (irma / "amended_consensus" / "Neg-KO-IAV_HA.fa").write_text(">x\nACGT\n")
    (irma / "amended_consensus" / "Neg-KO-IAV_NA.fa").write_text(">x\nACGT\n")
    names = tmp_path / "names.xlsx"
    assert process_sample_folder(irma, names, tmp_path) is None
    assert not names.exists()
    assert not (tmp_path / "GISAIDfasta").exists()

# CHANGE (2026-10-23): NEW. Naming a sample again replaces its rows in the
# run's GISAID files instead of adding a duplicate.
_SEQ = "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACC"
_SEGMENTS = ("HA", "NA", "PB2", "MP")


def _irma_sample(run_dir: Path, name: str, ha: str = "H1", na: str = "N1") -> Path:
    """Create ``<run_dir>/<name>/irma_results/<name>`` with an A/<ha><na> consensus."""
    irma = run_dir / name / "irma_results" / name
    (irma / "amended_consensus").mkdir(parents=True)
    (irma / f"A_HA_{ha}.fasta").write_text(f">A_HA_{ha}\n{_SEQ}\n")
    (irma / f"A_NA_{na}.fasta").write_text(f">A_NA_{na}\n{_SEQ}\n")
    for seg in _SEGMENTS:
        (irma / "amended_consensus" / f"{name}_{seg}.fa").write_text(f">{name}_{seg}\n{_SEQ}\n")
    return irma


def _run_tables(run_dir: Path):
    import pandas as pd
    out = run_dir / "GISAIDfasta"
    submit = pd.read_excel(out / "submit_info.xlsx")
    fastas = pd.read_excel(out / "fastas.xlsx", header=None)[0].tolist()
    return submit, [h for h in fastas if str(h).startswith(">")]


def test_rerun_replaces_rows_instead_of_duplicating(tmp_path: Path):
    pytest.importorskip("fastaparser")
    names = tmp_path / "names.xlsx"
    vi = _irma_sample(tmp_path, "2026-00418VI")
    other = _irma_sample(tmp_path, "26-057848", "H3", "N2")
    assert process_sample_folder(other, names, tmp_path) == "A/Finland/1/2026_H3N2"

    for _ in range(3):   # e.g. gisaid_naming_rerun.py --sample 2026-00418VI, run three times
        assert process_sample_folder(vi, names, tmp_path, rename_existing=True) == "A/Finland/2/2026_H1N1pdm09"

    submit, headers = _run_tables(tmp_path)
    assert sorted(submit["Isolate_name"]) == ["A/Finland/1/2026", "A/Finland/2/2026"]
    assert len(headers) == 2 * len(_SEGMENTS)
    assert len(set(headers)) == len(headers)


def test_force_new_number_removes_rows_under_the_old_name(tmp_path: Path):
    pytest.importorskip("fastaparser")
    names = tmp_path / "names.xlsx"
    vi = _irma_sample(tmp_path, "2026-00418VI")
    assert process_sample_folder(vi, names, tmp_path) == "A/Finland/1/2026_H1N1pdm09"
    assert process_sample_folder(vi, names, tmp_path, force_reassign=True) == "A/Finland/2/2026_H1N1pdm09"

    submit, headers = _run_tables(tmp_path)
    assert submit["Isolate_name"].tolist() == ["A/Finland/2/2026"]
    assert not any(h.startswith(">A/Finland/1/2026_") for h in headers)


# CHANGE (2026-10-23): NEW. A sample named in an earlier run is not named again.
def test_already_named_sample_is_skipped_in_a_later_run(tmp_path: Path):
    pytest.importorskip("fastaparser")
    import pandas as pd
    names = tmp_path / "names.xlsx"
    first_run, later_run = tmp_path / "RESULTS_20261014", tmp_path / "RESULTS_20261023"
    assert process_sample_folder(_irma_sample(first_run, "26-029774"), names, first_run) == \
        "A/Finland/1/2026_H1N1pdm09"
    before = pd.read_excel(names)

    resequenced = _irma_sample(later_run, "26-029774")
    assert process_sample_folder(resequenced, names, later_run) is None
    assert not (later_run / "GISAIDfasta").exists()          # nothing to submit again
    pd.testing.assert_frame_equal(pd.read_excel(names), before)

    # Asked for explicitly: same number, files in the later run.
    assert process_sample_folder(resequenced, names, later_run, rename_existing=True) == \
        "A/Finland/1/2026_H1N1pdm09"
    assert (later_run / "GISAIDfasta" / "submit_info.xlsx").exists()