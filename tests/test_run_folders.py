"""Results go to ``<output>/RESULTS_<YYYYMMDD>/``, dated when the pipeline starts."""
import gzip
import os
from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

import run_pipeline

FAKE_TOOLS = Path(__file__).parent / "fake_tools"
RUN = "261006_M01234_0042_000000000-ABCDE"


def _fastq(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("@r\nACGT\n+\nIIII\n")


class _Clock:
    """Stand-in for ``run_pipeline.datetime`` with a settable ``now()``."""

    current = datetime(2026, 10, 6, 23, 59)

    @classmethod
    def now(cls) -> datetime:
        return cls.current


@pytest.fixture
def pipeline_env(tmp_path, monkeypatch):
    pytest.importorskip("fastaparser")
    pytest.importorskip("openpyxl")
    monkeypatch.setenv("PATH", f"{FAKE_TOOLS}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_CONDA_LOG", str(tmp_path / "conda_calls.log"))
    monkeypatch.setattr(run_pipeline, "datetime", _Clock)
    (tmp_path / "config.json").write_text(
        '{"fastp_conda_env": "qc", "irma_conda_env": "irma_env", "qc_conda_env": "qc",'
        f' "logs_dir": "{tmp_path / "logs"}", "archive_dir": "{tmp_path / "archive"}",'
        ' "file_stable_wait_seconds": 0}'
    )
    return tmp_path


def _run(tmp_path: Path) -> None:
    run_pipeline.main(["-i", str(tmp_path / "input"), "-o", str(tmp_path / "output"),
                       "--config", str(tmp_path / "config.json"), "--max-workers", "1"])


def test_results_folder_is_named_after_the_start_date(pipeline_env):
    tmp = pipeline_env
    bc = tmp / "input" / "influenza" / "miseq" / RUN / "Data" / "Intensities" / "BaseCalls"
    for sample in ("26-040001_S1", "26-040002_S2"):
        for r in (1, 2):
            _fastq(bc / f"{sample}_L001_R{r}_001.fastq.gz")

    _Clock.current = datetime(2026, 10, 6, 23, 59)
    _run(tmp)

    assert [p.name for p in (tmp / "output").iterdir()] == ["RESULTS_20261006"]
    assert sorted(p.name for p in (tmp / "archive" / "RESULTS_20261006").iterdir()) == ["26-040001", "26-040002"]
    irma = pd.read_excel(tmp / "output" / "RESULTS_20261006" / "IRMA_analysis.xlsx")
    assert sorted(irma["Sample name"]) == ["26-040001", "26-040002"]
    assert (tmp / "logs" / "RESULTS_20261006" / "26-040001" / ".IRMA_DONE").exists()


def test_later_invocation_uses_its_own_date_and_log(pipeline_env):
    """Known trade-off of date-named folders: a sample completed by a later.
    
    invocation after midnight (e.g. its R2 arrived late) lands in the next
    day's folder. Each invocation also keeps its own log file.
    
    If the pipeline is run at 23:59 and again at 00:01, the first invocation
    will create RESULTS_20261006 and the second RESULTS_20261007. The first
    invocation will only see the R1 of 26-040002, and the second invocation
    will see the R2 of 26-040002. Each invocation will have its own log file.
    """
    tmp = pipeline_env
    bc = tmp / "input" / "influenza" / "miseq" / RUN / "Data" / "Intensities" / "BaseCalls"
    for r in (1, 2):
        _fastq(bc / f"26-040001_S1_L001_R{r}_001.fastq.gz")
    _fastq(bc / "26-040002_S2_L001_R1_001.fastq.gz")      # R2 still copying

    _Clock.current = datetime(2026, 10, 6, 23, 59)
    _run(tmp)
    _fastq(bc / "26-040002_S2_L001_R2_001.fastq.gz")      # arrives after midnight
    _Clock.current = datetime(2026, 10, 7, 0, 1)
    _run(tmp)

    assert sorted(p.name for p in (tmp / "output").iterdir()) == ["RESULTS_20261006", "RESULTS_20261007"]
    assert [p.name for p in (tmp / "archive" / "RESULTS_20261006").iterdir()] == ["26-040001"]
    assert [p.name for p in (tmp / "archive" / "RESULTS_20261007").iterdir()] == ["26-040002"]
    assert sorted(p.name for p in (tmp / "logs").glob("pipeline_*.log")) == \
        ["pipeline_20261006_235900.log", "pipeline_20261007_000100.log"]