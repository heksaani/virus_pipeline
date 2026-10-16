"""Generic per-sample workflow shared by every species and platform.
 
Every sample passes through the same three stages::
 
    platform read preparation  ->  IRMA assembly  ->  species post-IRMA step
 
Species-specific behaviour (IRMA module choice, post-IRMA processing) is
defined in ``species.py``. Platform-specific read preparation (Illumina vs
MinION) is defined in ``pre_irma.py``.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from irma import irma_external_config, run_irma
from pre_irma import check_env, run_minion_pre_irma_steps, run_pre_irma_steps

if TYPE_CHECKING:
    from species import SpeciesDef

logger = logging.getLogger(__name__)

@dataclass
class SampleRun:
    """Context for one sample, passed to the species' post-IRMA step.
 
    Attributes:
        sample: Sample name.
        platform: Sequencing platform, ``"illumina"`` or ``"minion"``.
        run_dir: Results folder for the run, ``<output>/<RUN_FOLDER>_RESULTS``.
        sample_dir: Results folder for this sample, ``run_dir/<sample>``.
        irma_dir: IRMA output folder, ``sample_dir/irma_results/<sample>``.
        reads: Read files given to IRMA: ``[R1, R2]`` for Illumina,
            ``[merged]`` for MinION.
        log_path: Per-sample log file.
        config: Pipeline configuration.
    """
    sample: str
    platform: str
    run_dir: Path
    sample_dir: Path
    irma_dir: Path
    reads: list[str]
    log_path: Path
    config: dict

    @property
    def amended_consensus(self) -> Path:
        """Folder holding IRMA's amended consensus sequences."""
        return self.irma_dir / "amended_consensus"


def prepare_reads(platform: str, sample: str, files_sorted: list, sample_dir: Path,
                  log_path: Path, config: dict) -> list[str] | None:
    """Run platform-specific QC and read preparation for one sample.
 
    Args:
        platform: ``"illumina"`` or ``"minion"``.
        sample: Sample name.
        files_sorted: Raw read files for the sample, in sorted order.
        sample_dir: Results folder for this sample.
        log_path: Per-sample log file.
        config: Pipeline configuration.
 
    Returns:
        Read file paths to give to IRMA, or ``None`` if a required step failed.
 
    Raises:
        ValueError: If ``platform`` is not ``"illumina"`` or ``"minion"``.
    """
    if platform == "illumina":
        return run_pre_irma_steps(sample, files_sorted, sample_dir, log_path, config)
    if platform == "minion":
        merged = run_minion_pre_irma_steps(sample, files_sorted, sample_dir, log_path, config)
        return [str(merged)] if merged else None
    raise ValueError(f"Unknown platform {platform!r}")


def run_sample_workflow(species: SpeciesDef, platform: str, sample: str, files_sorted: list,
                        run_dir: Path, log_path: Path, config: dict, cli_args=None) -> bool:
    """Run one sample through read preparation, IRMA and post-processing.
 
    A sample where IRMA ran but assembled nothing (e.g. a negative control)
    counts as finished: the species' post-IRMA step records the result and
    returns ``True``.
 
    Args:
        species: Species definition supplying the IRMA module and the
            post-IRMA step.
        platform: ``"illumina"`` or ``"minion"``.
        sample: Sample name.
        files_sorted: Raw read files for the sample, in sorted order.
        run_dir: Results folder for the run; the sample folder is created
            inside it.
        log_path: Per-sample log file.
        config: Pipeline configuration.
        cli_args: Not used. Kept so existing callers keep working.
 
    Returns:
        ``True`` if the sample completed, so ``run_pipeline.process_sample()``
        can archive its raw files and mark it done. ``False`` if read
        preparation, the IRMA environment check or IRMA itself failed.
    """
    logger.info("Starting %s %s workflow for %s", species.code, platform, sample)
    sample_dir = run_dir / sample
    sample_dir.mkdir(parents=True, exist_ok=True)

    reads = prepare_reads(platform, sample, files_sorted, sample_dir, log_path, config)
    if reads is None:
        return False

    irma_env = check_env("irma_conda_env", config)
    if not irma_env:
        return False
    irma_dir = sample_dir / "irma_results" / sample
    irma_dir.parent.mkdir(parents=True, exist_ok=True)
    if not run_irma(irma_env, reads, module=species.irma_module(platform, config),
                    sample_dir=irma_dir, log_path=log_path,
                    extra_conf=irma_external_config(config, species.code)):
        logger.error("IRMA failed for sample %s, skipping post-processing.", sample)
        return False
    logger.info("IRMA completed successfully for sample %s", sample)

    ctx = SampleRun(sample=sample, platform=platform, run_dir=run_dir, sample_dir=sample_dir,
                    irma_dir=irma_dir, reads=reads, log_path=log_path, config=config)
    return species.post_irma(ctx)
