"""End-to-end MinION runs (influenza and SARS-CoV-2)
on a MinKNOW-style run folder, with tests/fake_tools/conda standing in for
the real tools."""
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
FAKE_TOOLS = Path(__file__).parent / "fake_tools"
RUN = "20261006_1203_MN41234_FAX12345_abcd1234"


def _results_dir(tmp_path: Path) -> Path:
    """Return the run's single ``RESULTS_<YYYYMMDD>`` folder under ``output/``."""
    (results,) = (tmp_path / "output").iterdir()
    assert results.name.startswith("RESULTS_"), results.name
    return results


@pytest.fixture(autouse=True)
def _needs_python_tools():
    pytest.importorskip("fastaparser")
    pytest.importorskip("openpyxl")


def _reads(n: int, tag: str) -> str:
    return "".join(f"@{tag}{i}\nACGTACGT\n+\nIIIIIIII\n" for i in range(n))


def _minknow_run(tmp_path: Path, species: str) -> Path:
    """Two barcodes with pass + fail reads, gzipped and plain, plus unclassified."""
    run = tmp_path / "input" / species / "minion" / RUN
    for bc in ("barcode01", "barcode02"):
        for kind, n_files in (("pass", 2), ("fail", 1)):
            d = run / f"fastq_{kind}" / bc
            d.mkdir(parents=True)
            for i in range(n_files):
                with gzip.open(d / f"FAX12345_{kind}_{bc}_abcd1234_{i}.fastq.gz", "wt") as fh:
                    fh.write(_reads(3, f"{bc}_{kind}{i}_"))
        (run / "fastq_pass" / bc / "older_minknow_plain.fastq").write_text(_reads(2, f"{bc}_plain_"))
    unclassified = run / "fastq_pass" / "unclassified"
    unclassified.mkdir()
    for i in range(3):
        with gzip.open(unclassified / f"FAX12345_pass_unclassified_abcd1234_{i}.fastq.gz", "wt") as fh:
            fh.write(_reads(1, "u"))
    (run / "final_summary_FAX12345_abcd1234.txt").write_text("protocol_run_id=abcd1234\n")
    return run


def _config(tmp_path: Path) -> Path:
    (tmp_path / "nextclade_ds/flu_h3n2_ha").mkdir(parents=True)
    (tmp_path / "sars-cov-2-dataset").mkdir()
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps({
        "kraken_conda_env": "kraken", "kraken-db": str(tmp_path / "krakendb"),
        "nanoplot_conda_env": "nanoplot", "irma_conda_env": "irma_env",
        "snippy_conda_env": "snippy", "sars-ref": str(tmp_path / "ref.gbk"),
        "nextclade_conda_env": "nextclade", "nextclade_ref": str(tmp_path / "nextclade_ds"),
        "nextclade_ref_covid": str(tmp_path / "sars-cov-2-dataset"),
        "gisaid_naming_file": str(tmp_path / "names.xlsx"), "multiqc_conda_env": "qc",
        "logs_dir": str(tmp_path / "logs"), "archive_dir": str(tmp_path / "archive"),
        "file_stable_wait_seconds": 0,
    }))
    return cfg


def _run(tmp_path: Path) -> list[str]:
    env = dict(os.environ, PATH=f"{FAKE_TOOLS}{os.pathsep}{os.environ['PATH']}",
               FAKE_CONDA_LOG=str(tmp_path / "conda_calls.log"))
    proc = subprocess.run([sys.executable, str(REPO / "run_pipeline.py"), "-i", str(tmp_path / "input"),
                           "-o", str(tmp_path / "output"), "--config", str(_config(tmp_path)),
                           "--max-workers", "1"], env=env, capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr
    return (tmp_path / "conda_calls.log").read_text().splitlines()


def _count_reads(path: Path) -> int:
    with gzip.open(path, "rt") as fh:
        return sum(1 for line in fh if line.startswith("@"))


def test_minion_influenza_run(tmp_path: Path):
    run = _minknow_run(tmp_path, "influenza")
    calls = _run(tmp_path)
    run_dir = _results_dir(tmp_path)

    # one sample per barcode; unclassified skipped
    assert sorted(p.name for p in run_dir.iterdir() if p.name.startswith(("barcode", "unclassified"))) == \
        ["barcode01", "barcode02"]

    # pass + fail + plain reads merged into one gzipped file: 2*3 + 1*3 + 2 = 11
    merged = run_dir / "barcode01" / "barcode01.fastq.gz"
    assert _count_reads(merged) == 11

    irma_calls = [c for c in calls if " IRMA " in f" {c} "]
    assert len(irma_calls) == 2
    assert all(" IRMA FLU-minion " in c and c.count(".fastq.gz") == 1 for c in irma_calls)
    kraken_calls = [c for c in calls if "kraken2" in c]
    assert kraken_calls and not any("--paired" in c for c in kraken_calls)
    assert (run_dir / "barcode01" / "nanoplot" / "NanoStats.txt").exists()

    irma = pd.read_excel(run_dir / "IRMA_analysis.xlsx")
    assert sorted(irma["Sample name"]) == ["barcode01", "barcode02"]
    assert (run_dir / "segments" / "HA" / "barcode01_HA.fasta").exists()
    nextclade = pd.read_csv(run_dir / "nextclade_results.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert sorted(nextclade["sample_name"]) == ["barcode01", "barcode02"]

    # no GISAID naming for barcode names
    assert not (run_dir / "GISAIDfasta").exists()
    assert not (tmp_path / "names.xlsx").exists()

    # archived with pass/fail kept apart; unclassified left in place
    archived = tmp_path / "archive" / run_dir.name / "barcode01"
    assert len(list((archived / "fastq_pass").iterdir())) == 3
    assert len(list((archived / "fastq_fail").iterdir())) == 1
    assert not (run / "fastq_pass" / "barcode01").exists()
    assert len(list((run / "fastq_pass" / "unclassified").iterdir())) == 3
    assert (tmp_path / "logs" / run_dir.name / "barcode01" / ".IRMA_DONE").exists()


def test_minion_covid_run(tmp_path: Path):
    _minknow_run(tmp_path, "covid")
    calls = _run(tmp_path)
    run_dir = _results_dir(tmp_path)

    # IRMA has no "CoV-minion" config; the ONT config is minion-long-reads
    assert sum(" IRMA CoV-minion-long-reads " in f" {c} " for c in calls) == 2
    assert not any(" snippy " in f" {c} " for c in calls)   # short-read tool, not for ONT

    nextclade = pd.read_csv(run_dir / "nextclade_results_COV.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert sorted(nextclade["sample_name"]) == ["barcode01", "barcode02"]
    summary = pd.read_excel(run_dir / "sample_summary_report.xlsx").set_index("Sample name")
    assert summary.loc["barcode01", "Cov average coverage"] == 400
    assert (tmp_path / "logs" / run_dir.name / "barcode02" / ".IRMA_DONE").exists()