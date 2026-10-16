"""Run IRMA (Iterative Refinement Meta-Assembler) for one sample."""

import logging
import time
from pathlib import Path

from pre_irma import run_command_capture

logger = logging.getLogger(__name__)


def run_irma(conda_env: str,
             fastqs: list[str] | list[Path],
             module: str,
             sample_dir: Path,
             log_path: Path,
             dry_run: bool = False,
             extra_conf: Path | None = None) -> bool:
    """Assemble one sample with IRMA.

    Args:
        conda_env: Conda environment with IRMA.
        fastqs: The sample's reads: ``[R1, R2]`` (Illumina) or one merged
            file (MinION).
        module: IRMA module, e.g. ``"FLU"`` or ``"CoV-minion-long-reads"``.
        sample_dir: IRMA output folder for the sample.
        log_path: Log file to append IRMA's output to.
        dry_run: Log the command without running it.
        extra_conf: File passed to IRMA as ``--external-config``, if any.

    Returns:
        ``True`` if IRMA finished successfully (always ``True`` for a dry
        run), ``False`` if it failed.
    """
    sample_dir = Path(sample_dir)
    fastq_paths = [str(f) for f in fastqs]
    logger.info("IRMA start: module=%s sample_dir=%s fastqs=%s", module, sample_dir, fastq_paths)

    cmd = ["conda", "run", "-n", conda_env, "IRMA", module, *fastq_paths, str(sample_dir)]
    if extra_conf:
        cmd += ["--external-config", str(extra_conf)]

    if dry_run:
        logger.info("DRY RUN IRMA: %s", " ".join(cmd))
        return True

    logger.debug("IRMA command: %s", " ".join(cmd))
    start = time.monotonic()
    success = run_command_capture(cmd, log_path)
    elapsed = time.monotonic() - start

    if success:
        logger.info("IRMA OK: module=%s sample_dir=%s (%.1fs)", module, sample_dir, elapsed)
    else:
        logger.error("IRMA FAILED: module=%s sample_dir=%s (%.1fs), see %s",
                     module, sample_dir, elapsed, log_path)
    return success


def irma_external_config(config: dict, species_code: str | None = None) -> Path:
    """Find the file to pass to IRMA as ``--external-config`` for a species.

    The first of these that applies is used:

    1. ``config["irma_external_config_<code>"]``, e.g. ``irma_external_config_flu``.
    2. ``custom-config-<code>.sh`` next to this file, e.g. ``custom-config-flu.sh``.
    3. ``config["irma_external_config"]``.
    4. ``custom-config.sh`` next to this file.

    Per-species files keep species-only settings (such as influenza's
    ``SEG_NUMBERS``) out of other IRMA modules. IRMA reads the external
    config after the module's own ``init.sh`` and config, so its values win.

    Args:
        config: Pipeline configuration.
        species_code: Species code, e.g. ``"FLU"``; ``None`` skips steps 1-2.

    Returns:
        Path of the config file to use.
    """
    here = Path(__file__).parent
    if species_code:
        code = species_code.lower()
        if config.get(f"irma_external_config_{code}"):
            return Path(config[f"irma_external_config_{code}"])
        if (here / f"custom-config-{code}.sh").exists():
            return here / f"custom-config-{code}.sh"
    return Path(config.get("irma_external_config") or here / "custom-config.sh")