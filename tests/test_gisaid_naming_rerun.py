"""CHANGE (2026-10-06): NEW. MinION barcode samples are never GISAID-named by the re-run script."""
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