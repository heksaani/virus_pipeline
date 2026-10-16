"""Species registry: the one place to add a new organism to the pipeline.

Adding a species (details and an RSV example in README "Adding a species"):

  1. Add a ``SpeciesDef`` to ``SPECIES`` below.
  2. Give it a post-IRMA step: ``genome_post_irma`` (Nextclade on a single
     genome) for most viruses, or your own function taking a ``SampleRun``.
  3. Add its config keys (Nextclade dataset, ...) to ``config.json``.
  4. Copy a run into ``<input>/<folder>/<sequencer>/`` and run the pipeline.

Everything else is shared and needs no changes: the pre-IRMA steps
(Illumina: FastQC/fastp/Kraken2, MinION: merge/NanoPlot/Kraken2), IRMA,
archiving, done markers, the summary report and MultiQC.

IRMA module names must exist in ``IRMA_RES/modules/<MODULE>/config`` of your
IRMA install (https://github.com/CDCgov/irma). Use ``"MODULE"`` alone, or
``"MODULE-<config>"`` where ``config/<MODULE-config>.sh`` or
``config/<config>.sh`` exists; otherwise IRMA stops with
"Configuration ... not found".
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from covid_pipeline import covid_post_irma
from influenza_pipeline import flu_post_irma
from workflow import SampleRun, run_sample_workflow

# Sequencer folder name (input/<species>/<sequencer>/) -> platform.
SEQUENCER_PLATFORMS = {
    "nextseq": "illumina",
    "miseq": "illumina",
    "minion": "minion",
}

WorkflowHandler = Callable[..., bool]


@dataclass(frozen=True)
class SpeciesDef:
    """How the pipeline handles one organism.

    Attributes:
        code: Short code used in logs, config keys and routing, e.g. ``"FLU"``.
        folder: Input folder name under ``<input_dir>/``, e.g. ``"influenza"``.
        irma_modules: IRMA module per platform. A platform missing here is not
            supported for this species, and its files are left in place.
        post_irma: Species-specific steps after IRMA. Takes a ``SampleRun``
            and returns ``True`` when the sample's results are recorded.
        segmented: ``True`` for segmented genomes (influenza), which get
            per-segment summary columns.
        report_prefix: Summary-report column prefix for single-genome species,
            e.g. ``"Cov"`` gives ``"Cov N%"``, ``"Cov average coverage"``, ...
        aliases: Other folder names accepted for this species.
    """

    code: str
    folder: str
    irma_modules: dict[str, str]
    post_irma: Callable[[SampleRun], bool]
    segmented: bool = False
    report_prefix: str = ""
    aliases: tuple[str, ...] = field(default_factory=tuple)

    def irma_module(self, platform: str, config: dict) -> str:
        """Return the IRMA module to use for ``platform``.

        A lab can override the default with the config key
        ``irma_module_<code>_<platform>``, e.g. ``irma_module_flu_minion``.

        Args:
            platform: ``"illumina"`` or ``"minion"``.
            config: Pipeline configuration.

        Returns:
            IRMA module name, e.g. ``"FLU"`` or ``"CoV-minion-long-reads"``.

        Raises:
            ValueError: If the species does not support ``platform``. A config
                override cannot add a platform, because
                ``get_workflow_handler()`` only offers platforms listed in
                ``irma_modules``.
        """
        if platform not in self.irma_modules:
            raise ValueError(f"{self.code} has no IRMA module for platform {platform!r}")
        return config.get(f"irma_module_{self.code.lower()}_{platform}", self.irma_modules[platform])


SPECIES: list[SpeciesDef] = [
    SpeciesDef(
        code="FLU",
        folder="influenza",
        irma_modules={"illumina": "FLU", "minion": "FLU-minion"},
        post_irma=flu_post_irma,
        segmented=True,
    ),
    SpeciesDef(
        code="CoV",
        folder="covid",
        aliases=("sars-cov-2",),
        # IRMA has no "CoV-minion" config; its ONT configs are
        # minion-long-reads, minion-long-reads-notmp and minion-sensitive.
        irma_modules={"illumina": "CoV", "minion": "CoV-minion-long-reads"},
        post_irma=covid_post_irma,
        report_prefix="Cov",
    ),
]

def species_by_folder(folder: str) -> SpeciesDef | None:
    """Find a species by its input folder name or an alias (case-insensitive).

    Args:
        folder: Folder name under the input folder, e.g. ``"influenza"``.

    Returns:
        The matching ``SpeciesDef``, or ``None`` if no species matches.
    """
    folder = folder.lower()
    return next((s for s in SPECIES if folder == s.folder or folder in s.aliases), None)

def species_by_code(code: str | None) -> SpeciesDef | None:
    """Find a species by its exact code.

    Args:
        code: Species code, e.g. ``"FLU"``.

    Returns:
        The matching ``SpeciesDef``, or ``None`` if ``code`` is ``None`` or
        unknown.
    """
    return next((s for s in SPECIES if s.code == code), None)


def genome_report_prefixes() -> list[str]:
    """List the summary-report column prefixes of all single-genome species.

    Returns:
        Prefixes such as ``["Cov"]``, in registry order.
    """
    return [s.report_prefix for s in SPECIES if s.report_prefix]

def get_workflow_handler(species_code: str | None, seq_type: str | None) -> WorkflowHandler | None:
    """Return the workflow function for a species and sequencer pair.

    Args:
        species_code: Species code, e.g. ``"FLU"``.
        seq_type: Sequencer folder name, e.g. ``"nextseq"`` or ``"minion"``.

    Returns:
        A function with the signature ``run_pipeline.process_sample()``
        calls, ``handler(sample, files_sorted, run_dir, log_path, config,
        cli_args)``, returning ``True`` when the sample completed. ``None``
        if the species is unknown or does not support that sequencer.
    """
    species = species_by_code(species_code)
    platform = SEQUENCER_PLATFORMS.get(seq_type or "")
    if species is None or platform is None or platform not in species.irma_modules:
        return None

    def handler(sample: str, files_sorted: list, run_dir: Path, log_path: Path,
                config: dict, cli_args=None) -> bool:
        return run_sample_workflow(species, platform, sample, files_sorted,
                                   run_dir, log_path, config, cli_args)

    handler.__name__ = handler.__qualname__ = f"{species.code}_{platform}_workflow"
    return handler