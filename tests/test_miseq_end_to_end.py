"""End-to-end MiSeq smoke test.

Runs run_pipeline.py on a fake MiSeq run folder with `conda` replaced by
tests/fake_tools/conda, which writes small outputs shaped like real
fastqc / fastp / kraken2 / IRMA / nextclade output (IRMA amended consensus
uses the lab's SEG_NUMBERS names, or IRMA's stock numbering, plus
.pad.fa twins). This checks
the wiring between steps, not the bioinformatics.
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
FAKE_TOOLS = Path(__file__).parent / "fake_tools"


RUN_FOLDER = "261006_M01234_0042_000000000-ABCDE"


def _results_dir(tmp_path: Path) -> Path:
    """Return the run's single ``RESULTS_<YYYYMMDD>`` folder under ``output/``."""
    (results,) = (tmp_path / "output").iterdir()
    assert results.name.startswith("RESULTS_"), results.name
    return results


def _fastq(path: Path):
    import gzip
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("@r1\nACGTACGTAC\n+\nIIIIIIIIII\n")


def _setup(tmp_path: Path):
    bc = tmp_path / f"input/influenza/miseq/{RUN_FOLDER}/Data/Intensities/BaseCalls"
    for s in ("26-040001_S1", "Undetermined_S0"):
        for r in (1, 2):
            _fastq(bc / f"{s}_L001_R{r}_001.fastq.gz")
    (tmp_path / "nextclade_ds/flu_h3n2_ha").mkdir(parents=True)
    names = tmp_path / "names.xlsx"
    pd.DataFrame({"THL näytenimi": ["26-000001"], "juokseva numero": [139],
                  "GISAID nimi": ["A/Finland/139/2026_H3N2"], "huomoita": [None]}).to_excel(names, index=False)
    config = {
        "fastp_conda_env": "qc", "qc_conda_env": "qc", "kraken_conda_env": "kraken",
        "kraken-db": str(tmp_path / "krakendb"), "irma_conda_env": "irma_env",
        "nextclade_conda_env": "nextclade", "nextclade_ref": str(tmp_path / "nextclade_ds"),
        "gisaid_naming_file": str(names), "multiqc_conda_env": "qc",
        "logs_dir": str(tmp_path / "logs"), "archive_dir": str(tmp_path / "archive"),
        "file_stable_wait_seconds": 0,
    }
    (tmp_path / "config.json").write_text(json.dumps(config))
    return bc


def _run(tmp_path: Path, extra_env=None, *args):
    env = dict(os.environ, PATH=f"{FAKE_TOOLS}{os.pathsep}{os.environ['PATH']}",
               FAKE_CONDA_LOG=str(tmp_path / "conda_calls.log"), **(extra_env or {}))
    return subprocess.run(
        [sys.executable, str(REPO / "run_pipeline.py"), "-i", str(tmp_path / "input"),
         "-o", str(tmp_path / "output"), "--config", str(tmp_path / "config.json"),
         "--max-workers", "1", *args],
        env=env, capture_output=True, text=True, timeout=300,
    )


@pytest.fixture(autouse=True)
def _needs_python_tools():
    pytest.importorskip("fastaparser")
    pytest.importorskip("openpyxl")
    if shutil.which("bash") is None:
        pytest.skip("needs a POSIX shell")


@pytest.mark.parametrize("irma_naming", ["lab_seg_numbers", "numbers"])
def test_miseq_flu_run(tmp_path: Path, irma_naming):
    bc = _setup(tmp_path)
    proc = _run(tmp_path, {"FAKE_IRMA_NAMING": irma_naming})
    assert proc.returncode == 0, proc.stderr
    run_dir = _results_dir(tmp_path)

    # segments copied under their real names (not "1".."8" / "4.PAD")
    assert sorted(p.name for p in (run_dir / "segments").iterdir()) == \
        ["HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"]
    assert (run_dir / "segments/HA/26-040001_HA.fasta").exists()

    irma = pd.read_excel(run_dir / "IRMA_analysis.xlsx")
    row = irma.set_index("Sample name").loc["26-040001"]
    assert row["Subtype"] == "H3N2"
    assert row["HA N%"] == pytest.approx(0.99)  # from the .fa, not the .pad.fa

    names = pd.read_excel(tmp_path / "names.xlsx")
    assert names.iloc[-1]["GISAID nimi"] == "A/Finland/140/2026_H3N2"
    submit = pd.read_excel(run_dir / "GISAIDfasta/submit_info.xlsx")
    assert submit.loc[0, "Seq_Id (HA)"] == "A/Finland/140/2026_H3N2_HA"

    nextclade = pd.read_csv(run_dir / "nextclade_results.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert list(nextclade.columns[:3]) == ["sample_name", "Segmentti", "seqName"]
    assert nextclade["Segmentti"].tolist() == ["HA"]
    assert nextclade["substitutions"].tolist() == ["A12G,C34T;x"]  # commas/semicolons survive

    # sample archived and marked done; Undetermined left alone
    assert not (bc / "26-040001_S1_L001_R1_001.fastq.gz").exists()
    assert (tmp_path / "archive" / run_dir.name / "26-040001").is_dir()
    assert (tmp_path / "logs" / run_dir.name / "26-040001" / ".IRMA_DONE").exists()
    assert (bc / "Undetermined_S0_L001_R1_001.fastq.gz").exists()


def test_irma_failure_keeps_raw_files_for_retry(tmp_path: Path):
    bc = _setup(tmp_path)
    proc = _run(tmp_path, {"FAKE_IRMA_FAIL": "1"})
    assert proc.returncode == 0, proc.stderr
    run_dir_name = _results_dir(tmp_path).name
    assert (bc / "26-040001_S1_L001_R1_001.fastq.gz").exists()
    assert not (tmp_path / "logs" / run_dir_name / "26-040001" / ".IRMA_DONE").exists()


def test_dry_run_runs_no_tools(tmp_path: Path):
    bc = _setup(tmp_path)
    proc = _run(tmp_path, None, "--dry-run")
    assert proc.returncode == 0, proc.stderr
    calls = (tmp_path / "conda_calls.log")
    assert not calls.exists() or "run -n" not in calls.read_text()
    assert (bc / "26-040001_S1_L001_R1_001.fastq.gz").exists()