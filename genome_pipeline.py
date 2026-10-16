"""Post-IRMA step for single-genome (non-segmented) species.

SARS-CoV-2's Nextclade step generalised so
a new non-segmented virus (e.g. RSV, MERS-CoV) can reuse it from its
species.py entry without new code:

    post_irma=partial(genome_post_irma, nextclade_ref_key="nextclade_ref_rsv",
                      results_name="nextclade_results_RSV.tsv")

The per-sample summary row (Kraken, ambiguous sites, N%, coverage) is
written by run_pipeline for every species, so it is not repeated here.
"""
from __future__ import annotations

import logging
from pathlib import Path

from post_irma import run_nextclade_genome_analysis
from pre_irma import check_env
from workflow import SampleRun

logger = logging.getLogger(__name__)

def genome_post_irma(ctx: SampleRun, nextclade_ref_key: str, results_name: str,
                     summary_columns: list[str] | None = None) -> bool:
    """Nextclade on the consensus genome, if configured.

    Args:
        ctx: The sample (see workflow.SampleRun).
        nextclade_ref_key: config.json key holding the species' Nextclade dataset folder.
        results_name: Run-level TSV the sample's Nextclade row is appended to.
        summary_columns: Nextclade columns to keep (all if None).

    Returns:
        Always True: a missing env or dataset only logs a warning, and a
        sample with no consensus is a finished (negative) result.
    """
    nextclade_env = check_env("nextclade_conda_env", ctx.config)
    dataset = ctx.config.get(nextclade_ref_key)
    if not nextclade_env or not dataset:
        logger.warning("nextclade_conda_env or %s not set in config; skipping Nextclade for %s.",
                        nextclade_ref_key, ctx.sample)
        return True
    output_path = ctx.run_dir / "nextclade_results" / ctx.sample
    output_path.mkdir(parents=True, exist_ok=True)
    run_nextclade_genome_analysis(nextclade_env, Path(dataset), ctx.run_dir, ctx.amended_consensus,
                                  output_path, ctx.log_path, results_name=results_name,
                                  summary_columns=summary_columns, sample=ctx.sample)
    return True
