"""Influenza steps that run after IRMA.

Read preparation and IRMA itself are shared by all species (``workflow.py``).
After IRMA, an influenza sample goes through:

1. Copy each consensus segment to the run-level ``segments/<SEGMENT>/`` folder.
2. Add a row to ``IRMA_analysis.xlsx`` (N% per segment, read counts).
3. GISAID naming (Illumina only; MinION samples keep their barcode names).
4. Nextclade per segment.

``flu_post_irma`` is the entry point registered in ``species.py``.
"""
import logging
from pathlib import Path

import fastaparser

from gisaid_naming import process_sample_folder as run_gisaid_naming
from post_irma import analyze_and_append_sample, run_nexclade_analysis
from pre_irma import check_env
from segments import FLU_SEGMENTS, consensus_segments
from workflow import SampleRun

logger = logging.getLogger(__name__)


def flu_post_irma(ctx: SampleRun) -> bool:
    """Run the influenza post-IRMA steps for one sample (entry point for ``species.py``).

    MinION samples are never GISAID-named, because barcode names are not THL
    sample names.

    Args:
        ctx: The sample's run context from ``workflow.run_sample_workflow``.

    Returns:
        ``True`` once the sample's results are recorded (see
        ``run_flu_post_irma``).
    """
    return run_flu_post_irma(ctx.sample, ctx.irma_dir.parent, ctx.run_dir, ctx.log_path, ctx.config,
                             gisaid_naming=ctx.platform == "illumina")


def write_segment_fastas(sample: str, irma_sample_dir: Path, run_dir: Path) -> int:
    """Copy each consensus segment to ``<run_dir>/segments/<SEGMENT>/<sample>_<SEGMENT>.fasta``.

    Only real consensus files are copied (not IRMA's N-padded ``*.pad.fa``),
    and both of IRMA's naming styles are understood (``<sample>_HA.fa`` and
    ``<sample>_4.fa``). Segments with an empty sequence are skipped.

    Args:
        sample: Sample name, used in the output file name and FASTA header.
        irma_sample_dir: IRMA output folder of the sample.
        run_dir: Results folder of the run.

    Returns:
        Number of segment files written.
    """
    written = 0
    for fasta, seg in consensus_segments(irma_sample_dir / "amended_consensus"):
        sequence = ""
        with open(fasta) as fasta_file:
            for record in fastaparser.Reader(fasta_file, parse_method="quick"):
                sequence = record.sequence

        if not sequence:
            logger.warning("Sample %s: %s has no consensus sequence for segment %s, skipping it.",
                           sample, fasta.name, seg)
            continue

        target = run_dir / "segments" / seg / f"{sample}_{seg}.fasta"
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(target, "w") as outfh:
            fastaparser.Writer(outfh).writefasta((f"{sample}_{seg}", sequence))
        logger.info("Wrote %s", target)
        written += 1
    return written


def run_flu_post_irma(sample: str, irma_results_dir: Path, run_dir: Path,
                      log_path: Path, config: dict, gisaid_naming: bool = True) -> bool:
    """Run segments, ``IRMA_analysis.xlsx``, GISAID naming and Nextclade for one sample.

    A sample where IRMA ran but assembled nothing (e.g. a negative control)
    is a finished result, not a failure: it gets a "No consensus" row in
    ``IRMA_analysis.xlsx``. Optional steps whose environment or config is
    missing are skipped with a warning. GISAID naming and Nextclade are
    independent, so a naming error does not stop Nextclade.

    Args:
        sample: Sample name.
        irma_results_dir: IRMA results folder that contains ``<sample>/``.
        run_dir: Results folder of the run.
        log_path: Log file to append tool output to.
        config: Pipeline configuration (``gisaid_naming_file``,
            ``nextclade_conda_env``, ``nextclade_ref``).
        gisaid_naming: ``False`` skips GISAID naming (MinION barcodes).

    Returns:
        ``True`` once the sample's results are recorded, including for a
        sample without a consensus.
    """
    for segment in FLU_SEGMENTS:
        (run_dir / "segments" / segment).mkdir(parents=True, exist_ok=True)

    irma_sample_dir = irma_results_dir / sample
    amended = irma_sample_dir / "amended_consensus"
    if not amended.exists():
        logger.warning("Sample %s has no amended_consensus folder; recording it as 'No consensus'.", sample)
        analyze_and_append_sample(irma_results_dir, sample, run_dir / "IRMA_analysis.xlsx")
        return True

    write_segment_fastas(sample, irma_sample_dir, run_dir)
    analyze_and_append_sample(irma_results_dir, sample, run_dir / "IRMA_analysis.xlsx")

    gisaid_excel = config.get("gisaid_naming_file")
    if not gisaid_naming:
        logger.info("GISAID naming not done for %s (MinION barcode name).", sample)
    elif not gisaid_excel:
        logger.warning("gisaid_naming_file not specified in config. Skipping GISAID naming step.")
    else:
        try:
            new_name = run_gisaid_naming(irma_sample_dir, Path(gisaid_excel), run_dir)
        except RuntimeError:
            logger.exception("GISAID naming failed for sample %s", sample)
        else:
            if new_name:
                logger.info("GISAID name for sample %s: %s", sample, new_name)
            else:
                logger.warning("Sample %s was not GISAID-named (not a THL sample name, or missing HA/NA).",
                               sample)

    nextclade_env = check_env("nextclade_conda_env", config)
    nextclade_ref = config.get("nextclade_ref")
    if not nextclade_env or not nextclade_ref:
        logger.warning("nextclade_conda_env or nextclade_ref not specified in config. Skipping Nextclade step.")
        return True
    output_path = run_dir / "nextclade_results" / sample
    output_path.mkdir(parents=True, exist_ok=True)
    run_nexclade_analysis(nextclade_env, Path(nextclade_ref), run_dir, amended, output_path, log_path,
                          sample_root=irma_sample_dir)
    return True