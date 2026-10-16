"""CHANGE (2026-10-06): NEW. process_sample() archiving / done-marker rules,
with the workflow handler replaced by a stub so no tools are needed."""
from argparse import Namespace
from pathlib import Path

import pytest

import run_pipeline


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
