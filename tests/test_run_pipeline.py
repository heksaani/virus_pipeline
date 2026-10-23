"""Script tests for run_pipeline.py, which is used to run the IRMA pipeline on a set of input files.
"""
from argparse import Namespace
from pathlib import Path

import pytest

import run_pipeline
from pre_irma import group_by_sample


def _make_sample(input_dir: Path, species="influenza", seq="miseq", name="26-040001_S1"):
    d = input_dir / species / seq / "run"
    d.mkdir(parents=True, exist_ok=True)
    files = []
    for r in (1, 2):
        f = d / f"{name}_L001_R{r}_001.fastq.gz"
        f.write_bytes(b"data")
        files.append(f)
    return files


@pytest.fixture
def env(tmp_path, monkeypatch):
    input_dir, run_dir = tmp_path / "in", tmp_path / "out" / "20261006_RESULTS"
    run_dir.mkdir(parents=True)
    config = {"archive_dir": str(tmp_path / "archive"), "file_stable_wait_seconds": 0}
    args = Namespace(logs_dir=tmp_path / "logs", force=False, dry_run=False)
    # Summary row needs no real outputs, but keep tests independent of it.
    monkeypatch.setattr(run_pipeline, "write_summary_row", lambda *a, **k: None)
    return input_dir, run_dir, config, args, tmp_path


def _stub_handler(monkeypatch, result):
    calls = []

    def handler(sample, files, run_dir, log, config, args):
        calls.append(sample)
        return result
    handler.__name__ = "stub"
    monkeypatch.setattr(run_pipeline, "get_workflow_handler", lambda sp, st: handler)
    return calls


def test_success_archives_and_marks_done(env, monkeypatch):
    input_dir, run_dir, config, args, tmp = env
    files = _make_sample(input_dir)
    _stub_handler(monkeypatch, True)
    run_pipeline.process_sample("26-040001", files, run_dir, input_dir, config, args)
    assert not any(f.exists() for f in files)
    assert (tmp / "archive" / run_dir.name / "26-040001" / files[0].name).exists()
    assert (tmp / "logs" / run_dir.name / "26-040001" / ".IRMA_DONE").exists()


def test_failed_workflow_leaves_files_and_no_marker(env, monkeypatch):
    input_dir, run_dir, config, args, tmp = env
    files = _make_sample(input_dir)
    _stub_handler(monkeypatch, False)
    run_pipeline.process_sample("26-040001", files, run_dir, input_dir, config, args)
    assert all(f.exists() for f in files)
    assert not (tmp / "logs" / run_dir.name / "26-040001" / ".IRMA_DONE").exists()


def test_unsupported_combination_is_left_in_place(env):
    input_dir, run_dir, config, args, _ = env
    files = _make_sample(input_dir, species="covid", seq="iseq")  # unknown sequencer
    run_pipeline.process_sample("26-040001", files, run_dir, input_dir, config, args)
    assert all(f.exists() for f in files)


def test_dry_run_changes_nothing(env, monkeypatch):
    input_dir, run_dir, config, args, _tmp = env
    args.dry_run = True
    files = _make_sample(input_dir)
    calls = _stub_handler(monkeypatch, True)
    run_pipeline.process_sample("26-040001", files, run_dir, input_dir, config, args)
    assert calls == []
    assert all(f.exists() for f in files)


def test_undetermined_is_skipped(env, monkeypatch):
    input_dir, run_dir, config, args, _ = env
    files = _make_sample(input_dir, name="Undetermined_S0")
    calls = _stub_handler(monkeypatch, True)
    run_pipeline.process_sample("Undetermined", files, run_dir, input_dir, config, args)
    assert calls == []
    assert all(f.exists() for f in files)


def test_basespace_download_folder_is_archived_and_removed(env, monkeypatch):
    """BaseSpace puts each sample in ``<sample>_L001-ds.<id>/``.

    The sample name comes from the file name, not the folder, and the
    emptied download folder is removed; the sequencer folder stays.
    """
    input_dir, run_dir, config, args, tmp = env
    miseq = input_dir / "influenza" / "miseq"
    download = miseq / "2026-00418VI_L001-ds.62d11019a0124daf989a947404652f18"
    download.mkdir(parents=True)
    files = []
    for r in (1, 2):
        f = download / f"2026-00418VI_S1_L001_R{r}_001.fastq.gz"
        f.write_bytes(b"data")
        files.append(f)

    assert list(group_by_sample(files)) == ["2026-00418VI"]
    _stub_handler(monkeypatch, True)
    assert run_pipeline.process_sample("2026-00418VI", files, run_dir, input_dir, config, args)

    assert (tmp / "archive" / run_dir.name / "2026-00418VI" / files[0].name).exists()
    assert not download.exists()
    assert miseq.is_dir()


@pytest.mark.parametrize("folder, sample, expected", [
    ("26-040001", "26-040001", True),                   # plain sample folder
    ("26-040001_L001-ds.abc123", "26-040001", True),    # BaseSpace download folder
    ("barcode01", "barcode01", True),                   # MinION barcode
    ("26-040001_L001-ds.abc123", "26-04000", False),    # other sample's prefix
    ("BaseCalls", "26-040001", False),                  # shared run folder
])
def test_is_sample_folder(folder, sample, expected):
    assert run_pipeline.is_sample_folder(Path(folder), sample) is expected