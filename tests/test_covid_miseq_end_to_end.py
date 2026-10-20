"""End-to-end MiSeq COVID run with fake tools
(see tests/fake_tools/conda): fastp -> Kraken2 -> IRMA CoV -> Snippy -> Nextclade.
"""
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
RUN = "261006_M01234_0043_000000000-COVID"


def _results_dir(tmp_path: Path) -> Path:
    """Return the run's single ``RESULTS_<YYYYMMDD>`` folder under ``output/``."""
    (results,) = (tmp_path / "output").iterdir()
    assert results.name.startswith("RESULTS_"), results.name
    return results


@pytest.fixture(autouse=True)
def _needs_python_tools():
    pytest.importorskip("fastaparser")
    pytest.importorskip("openpyxl")


def test_miseq_covid_run(tmp_path: Path):
    bc = tmp_path / f"input/covid/miseq/{RUN}/Data/Intensities/BaseCalls"
    bc.mkdir(parents=True)
    for r in (1, 2):
        with gzip.open(bc / f"26-050001_S1_L001_R{r}_001.fastq.gz", "wt") as fh:
            fh.write("@r\nACGT\n+\nIIII\n")
    (tmp_path / "sars-cov-2-dataset").mkdir()
    (tmp_path / "MN908947.gbk").write_text("LOCUS")
    (tmp_path / "config.json").write_text(json.dumps({
        "fastp_conda_env": "qc", "qc_conda_env": "qc", "kraken_conda_env": "kraken",
        "kraken-db": str(tmp_path / "krakendb"), "irma_conda_env": "irma_env",
        "snippy_conda_env": "snippy", "sars-ref": str(tmp_path / "MN908947.gbk"),
        "nextclade_conda_env": "nextclade", "nextclade_ref_covid": str(tmp_path / "sars-cov-2-dataset"),
        "multiqc_conda_env": "qc", "logs_dir": str(tmp_path / "logs"),
        "archive_dir": str(tmp_path / "archive"), "file_stable_wait_seconds": 0,
    }))
    env = dict(os.environ, PATH=f"{FAKE_TOOLS}{os.pathsep}{os.environ['PATH']}",
               FAKE_CONDA_LOG=str(tmp_path / "conda_calls.log"))
    proc = subprocess.run([sys.executable, str(REPO / "run_pipeline.py"), "-i", str(tmp_path / "input"),
                           "-o", str(tmp_path / "output"), "--config", str(tmp_path / "config.json"),
                                                      "--max-workers", "1"], check=False, env=env, capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr
    run_dir = _results_dir(tmp_path)

    calls = (tmp_path / "conda_calls.log").read_text()
    assert "IRMA CoV" in calls
    assert "snippy" in calls
    assert (run_dir / "26-050001" / "snippy" / "snps.tab").exists()

    nextclade = pd.read_csv(run_dir / "nextclade_results_COV.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert nextclade["sample_name"].tolist() == ["26-050001"]

    summary = pd.read_excel(run_dir / "sample_summary_report.xlsx").set_index("Sample name")
    row = summary.loc["26-050001"]
    assert row["Cov ambiguous count"] == 1          # one consensus site below 0.95
    assert row["Cov average coverage"] == 400
    assert row["Cov N%"] == pytest.approx(0.25)     # from the .fa, not the N-padded .pad.fa

    assert not any(bc.iterdir())                    # raw files archived
    assert (tmp_path / "logs" / run_dir.name / "26-050001" / ".IRMA_DONE").exists()