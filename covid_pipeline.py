"""SARS-CoV-2 post-IRMA steps.

Post-IRMA steps for SARS-CoV-2 include running Snippy (Illumina only) and Nextclade.

"""
import logging
from pathlib import Path

from genome_pipeline import genome_post_irma
from post_irma import NEXTCLADE_SUMMARY_COLUMNS_COVID
from pre_irma import check_env, run_snippy
from workflow import SampleRun

logger = logging.getLogger(__name__)


def covid_post_irma(ctx: SampleRun) -> bool:
    """Run post-IRMA steps for SARS-CoV-2.

    Args:
        ctx: SampleRun context.

    Returns:
        ``True`` if the sample was processed successfully, ``False`` if it was skipped (e.g., no consensus genome).
    """
    if ctx.platform == "illumina":
        snippy_env = check_env("snippy_conda_env", ctx.config)
        sars_ref = ctx.config.get("sars-ref")
        if snippy_env and sars_ref:
            r1, r2 = ctx.reads
            run_snippy(snippy_env, r1, r2, ref=Path(sars_ref),
                       out_dir=ctx.sample_dir / "snippy", log_path=ctx.log_path,
                       threads=ctx.config.get("snippy_threads", 4))
        else:
            logger.warning("snippy_conda_env or sars-ref missing; skipping Snippy for %s.", ctx.sample)

    return genome_post_irma(ctx, nextclade_ref_key="nextclade_ref_covid",
                            results_name="nextclade_results_COV.tsv",
                            summary_columns=NEXTCLADE_SUMMARY_COLUMNS_COVID)
