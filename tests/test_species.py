"""The species registry, and adding a species
without touching any workflow code.
"""
import gzip
import json
import os
from functools import partial
from pathlib import Path

import pandas as pd
import pytest

import irma
import run_pipeline
import species
from genome_pipeline import genome_post_irma

FAKE_TOOLS = Path(__file__).parent / "fake_tools"


def test_routing_comes_from_registry():
    assert species.get_workflow_handler("FLU", "miseq").__name__ == "FLU_illumina_workflow"
    assert species.get_workflow_handler("CoV", "minion").__name__ == "CoV_minion_workflow"
    assert species.get_workflow_handler("FLU", "iseq") is None      # unknown sequencer
    assert species.get_workflow_handler("XYZ", "miseq") is None     # unknown species


def test_irma_module_can_be_overridden_in_config():
    cov = species.species_by_code("CoV")
    assert cov.irma_module("minion", {}) == "CoV-minion-long-reads"
    assert cov.irma_module("minion", {"irma_module_cov_minion": "CoV-minion-sensitive"}) == "CoV-minion-sensitive"


def test_irma_external_config_per_species(tmp_path):
    here = Path(irma.__file__).parent
    assert irma.irma_external_config({}, "FLU") == here / "custom-config-flu.sh"   # has SEG_NUMBERS
    assert irma.irma_external_config({}, "CoV") == here / "custom-config.sh"
    assert irma.irma_external_config({"irma_external_config_cov": "/x/cov.sh"}, "CoV") == Path("/x/cov.sh")


def test_flu_external_config_keeps_lab_seg_numbers():
    text = (Path(irma.__file__).parent / "custom-config-flu.sh").read_text()
    assert 'SEG_NUMBERS="B_PB1:B_PB1,B_PB2:B_PB2,A_PB2:A_PB2,A_PB1:A_PB1,PA:PA,HA:HA,NP:NP,NA:NA,MP:MP,NS:NS"' in text


def test_new_species_needs_only_a_registry_entry(tmp_path, monkeypatch):
    """Register an RSV-like species and run a MiSeq sample through the real
    pipeline (fake tools): routing, IRMA module, Nextclade, summary columns,
    archiving and done marker all come from the registry entry.
    """
    pytest.importorskip("fastaparser")
    pytest.importorskip("openpyxl")
    rsv = species.SpeciesDef(
        code="RSV", folder="rsv",
        irma_modules={"illumina": "RSV", "minion": "RSV-minion"},
        post_irma=partial(genome_post_irma, nextclade_ref_key="nextclade_ref_rsv",
                          results_name="nextclade_results_RSV.tsv"),
        report_prefix="RSV",
    )
    monkeypatch.setattr(species, "SPECIES", species.SPECIES + [rsv])

    run = "261006_M01234_0050_000000000-RSVAA"
    bc = tmp_path / "input" / "rsv" / "miseq" / run
    bc.mkdir(parents=True)
    for r in (1, 2):
        with gzip.open(bc / f"26-060001_S1_L001_R{r}_001.fastq.gz", "wt") as fh:
            fh.write("@r\nACGT\n+\nIIII\n")
    (tmp_path / "rsv_dataset").mkdir()
    (tmp_path / "config.json").write_text(json.dumps({
        "fastp_conda_env": "qc", "irma_conda_env": "irma_env", "nextclade_conda_env": "nextclade",
        "nextclade_ref_rsv": str(tmp_path / "rsv_dataset"),
        "logs_dir": str(tmp_path / "logs"), "archive_dir": str(tmp_path / "archive"),
        "file_stable_wait_seconds": 0,
    }))
    monkeypatch.setenv("PATH", f"{FAKE_TOOLS}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_CONDA_LOG", str(tmp_path / "calls.log"))

    run_pipeline.main(["-i", str(tmp_path / "input"), "-o", str(tmp_path / "output"),
                       "--config", str(tmp_path / "config.json"), "--max-workers", "1"])

    (run_dir,) = (tmp_path / "output").iterdir()
    assert " IRMA RSV " in (tmp_path / "calls.log").read_text()
    nextclade = pd.read_csv(run_dir / "nextclade_results_RSV.tsv", sep="\t", dtype=str)
    assert nextclade["sample_name"].tolist() == ["26-060001"]
    summary = pd.read_excel(run_dir / "sample_summary_report.xlsx").set_index("Sample name")
    assert summary.loc["26-060001", "RSV average coverage"] == 400
    assert (tmp_path / "logs" / run_dir.name / "26-060001" / ".IRMA_DONE").exists()
    assert not any(bc.iterdir())   # archived