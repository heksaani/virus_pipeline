r"""Run the IRMA pipeline on every sample found in an input folder.

For each sample: skip checks -> species workflow (read preparation, IRMA,
species post-processing) -> summary-report row -> archive raw FASTQs ->
``.IRMA_DONE`` marker. Samples run in parallel. All results of one pipeline
run go to ``<output>/RESULTS_<YYYYMMDD>/`` (the date the pipeline started),
followed by one MultiQC report if any sample completed.

Usage::

    python run_pipeline.py -i /path/to/input -o /path/to/output \
        --config config.json [--logs-dir DIR] [--max-workers N] [--force] [--dry-run]

Expected input layout::

    <input>/<species folder>/<sequencer folder>/.../*.fastq.gz

Species folders are defined in species.py; sequencer folders are
``nextseq``, ``miseq`` and ``minion``.
"""
import argparse
import json
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

from post_irma import append_sample_row, report_header, run_multiqc
from pre_irma import (
    check_env,
    find_fastq_files,
    get_species_and_seq_type,
    group_by_sample,
    has_r1_r2,
    unstable_files,
)
from species import genome_report_prefixes, get_workflow_handler, species_by_code

logger = logging.getLogger(__name__)

# Sized for an 8-core VM: each tool call uses --threads 4 by default, so two
# samples at a time use all cores without oversubscribing. Override with
# --max-workers or config "max_parallel_samples".
DEFAULT_MAX_WORKERS = 2

DEFAULT_ARCHIVE_DIR = "/mnt/volume/storage"
DEFAULT_LOGS_DIR = "/mnt/volume/logs"

# Not real samples: bcl2fastq/BCL Convert's reads with no matching index, and
# MinKNOW's reads with no barcode. Skipped and left in place.
SKIP_SAMPLE_NAMES = {"Undetermined", "unclassified"}

def is_sample_folder(folder: Path, sample: str) -> bool:
    """Check whether an input folder belongs to one sample only.

    True for a folder named after the sample (``26-040001/``, MinION
    ``barcode01/``) or starting with ``<sample>_`` (a BaseSpace download
    folder such as ``26-040001_L001-ds.<id>/``). The ``_`` keeps
    ``26-04000`` from matching ``26-040001_...``.

    Args:
        folder: Folder that held the sample's raw FASTQs.
        sample: Sample name.

    Returns:
        ``True`` if the folder may be removed once it is empty.
    """
    return folder.name == sample or folder.name.startswith(f"{sample}_")


def archive_sample_files(files_sorted: list[str], sample: str, archive_root: Path) -> None:
    """Move a sample's raw FASTQs to ``<archive_root>/<sample>/``.

    MinION files go to ``<sample>/<fastq_pass|fastq_fail>/<file>``, because
    the pass and fail folders can hold files with the same name. Files that
    already exist at the destination are left in place with a warning. An
    input folder of the sample (see ``is_sample_folder``) is removed once it
    is empty.

    Args:
        files_sorted: The sample's raw FASTQ files.
        sample: Sample name.
        archive_root: Archive folder for this run.
    """
    sample_storage_dir = archive_root / sample
    sample_storage_dir.mkdir(parents=True, exist_ok=True)
    for f in files_sorted:
        src = Path(f)
        if src.parent.name == sample and src.parent.name.lower().startswith("barcode"):
            dest = sample_storage_dir / src.parent.parent.name / src.name
            dest.parent.mkdir(parents=True, exist_ok=True)
        else:
            dest = sample_storage_dir / src.name
        if dest.exists():
            logger.warning("Destination file already exists, skipping move: %s", dest)
            continue
        logger.info("Moving %s to %s", src, dest)
        source_parent = src.parent
        src.rename(dest)
        
        if is_sample_folder(source_parent, sample) and source_parent.is_dir() and not any(source_parent.iterdir()):
            logger.info("Removing empty sample folder %s", source_parent)
            source_parent.rmdir()


def write_summary_row(sample: str, sample_species: str | None, run_dir: Path) -> None:
    """Add or update the sample's row in ``run_dir/sample_summary_report.xlsx``.

    The columns a sample fills come from its species.py entry: per-segment
    columns for segmented species (influenza), ``"<prefix> ..."`` columns
    for single-genome species. Missing result files leave their columns empty.

    Args:
        sample: Sample name.
        sample_species: Species code, or ``None`` if unknown.
        run_dir: Results folder for the run.
    """
    species = species_by_code(sample_species)
    sample_root = run_dir / sample
    irma_dir = sample_root / "irma_results" / sample

    def existing(p: Path) -> Path | None:
        return p if p.exists() else None

    append_sample_row(
        sample,
        run_dir / "sample_summary_report.xlsx",
        segmented=bool(species and species.segmented),
        genome_prefix=species.report_prefix if species else None,
        kraken_report_path=existing(sample_root / "kraken" / f"{sample}_kraken_report.txt"),
        irma_tables_dir=existing(irma_dir / "tables"),
        amended_consensus_dir=existing(irma_dir / "amended_consensus"),
        columns=report_header(genome_report_prefixes()),
    )

def process_sample(sample: str,
                   files: list[Path],
                   run_dir: Path,
                   input_dir: Path,
                   config: dict,
                   args: argparse.Namespace) -> bool:
    """Run every step for one sample.

    Order: skip checks -> workflow -> summary row -> archive raw files ->
    write ``.IRMA_DONE``. The sample is skipped, with its raw files left in
    ``input_dir`` for a later run, if it:

    * is not a real sample (``SKIP_SAMPLE_NAMES``),
    * is already marked done (unless ``--force``),
    * is an Illumina sample still missing R1 or R2,
    * has files that are still being written,
    * has no workflow for its species and sequencer,
    * is a dry run, or
    * its workflow failed or raised.

    Expected per-sample problems are logged, never raised, so one sample
    cannot stop the rest of the run.

    Args:
        sample: Sample name.
        files: The sample's raw FASTQ files.
        run_dir: Results folder of this pipeline run.
        input_dir: Pipeline input folder.
        config: Pipeline configuration.
        args: Parsed command-line arguments (``logs_dir``, ``force``, ``dry_run``).

    Returns:
        ``True`` if the sample completed in this call. ``main()`` runs
        MultiQC only if at least one sample completed.
    """
    if sample in SKIP_SAMPLE_NAMES:
        logger.info("Skipping %s (not a sample; left in place).", sample)
        return False

    sample_log_dir = Path(args.logs_dir).resolve() / run_dir.name / sample
    sample_log_dir.mkdir(parents=True, exist_ok=True)
    done_marker = sample_log_dir / ".IRMA_DONE"
    if done_marker.exists() and not args.force:
        logger.info("Sample %s already marked done (%s), skipping", sample, done_marker)
        return False

    files_sorted = sorted(map(str, files))
    sample_log = sample_log_dir / "pipeline.log"

    sample_species, sample_seq_type = get_species_and_seq_type(files_sorted, input_dir)

    # MinION barcodes are single-end, so only Illumina samples need a pair.
    if sample_seq_type != "minion":
        has_r1, has_r2 = has_r1_r2(files_sorted)
        if not (has_r1 and has_r2):
            logger.warning(
                "Sample %s does not yet have both R1 and R2 (R1=%s, R2=%s); "
                "leaving files in place for a later run instead of archiving them.",
                sample, has_r1, has_r2,
            )
            return False

    wait = float(config.get("file_stable_wait_seconds", 3.0))
    unstable = unstable_files([Path(f) for f in files_sorted], wait_seconds=wait)
    if unstable:
        logger.warning(
            "Sample %s has %d file(s) still being written, skipping this run: %s",
            sample, len(unstable), [str(u) for u in unstable[:5]],
        )
        return False

    handler = get_workflow_handler(sample_species, sample_seq_type)
    if handler is None:
        logger.warning("Not supported pipeline: %s %s for sample %s; leaving raw files in place.",
                       sample_species, sample_seq_type, sample)
        return False

    if args.dry_run:
        logger.info("DRY RUN: would run %s for sample %s on %s", handler.__name__, sample, files_sorted)
        return False

    logger.info("Sample %s: running workflow handler for species %s and sequencer type %s",
                sample, sample_species, sample_seq_type)
    try:
        ok = handler(sample, files_sorted, run_dir, sample_log, config, args)
    except Exception:
        logger.exception(
            "Unhandled error while processing sample %s; skipping this sample "
            "and continuing with the rest of the run. Raw files left in place "
            "in %s for inspection/retry.",
            sample, input_dir,
        )
        return False
    if not ok:
        logger.warning("Workflow for sample %s did not complete; raw files left in place for a retry. "
                       "See %s", sample, sample_log)
        return False

    # The summary row is a report, not a reason to keep re-running the
    # sample, so a failure here is logged and the sample still completes.
    try:
        write_summary_row(sample, sample_species, run_dir)
    except Exception:
        logger.exception("Could not add sample %s to sample_summary_report.xlsx", sample)

    archive_root = Path(config.get("archive_dir", DEFAULT_ARCHIVE_DIR)) / run_dir.name
    archive_sample_files(files_sorted, sample, archive_root)
    done_marker.touch()
    logger.info("Sample %s done.", sample)
    return True


def setup_logging(log_file: Path, with_thread_names: bool) -> None:
    """Send log messages to ``log_file`` (overwritten) and the console.

    Args:
        log_file: Run log file.
        with_thread_names: Add the worker thread name to each line; useful
            when samples run in parallel.
    """
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
    fmt_str = ("%(asctime)s - %(threadName)s - %(levelname)s - %(message)s" if with_thread_names
               else "%(asctime)s - %(levelname)s - %(message)s")
    fmt = logging.Formatter(fmt_str)
    for handler in (logging.FileHandler(log_file, mode="w"), logging.StreamHandler()):
        handler.setFormatter(fmt)
        root.addHandler(handler)
    root.setLevel(logging.INFO)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the command-line arguments.

    Args:
        argv: Arguments to parse; ``None`` uses ``sys.argv``.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Run IRMA across samples and postprocess results")
    parser.add_argument("-i", "--input", required=True, type=Path,
                        help="Input folder: <species>/<sequencer>/<run folder>/...")
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="Output folder; results go to <output>/RESULTS_<YYYYMMDD>/")
    parser.add_argument("--logs-dir", type=Path, default=None,
                        help="Log folder (default: config 'logs_dir', "
                             f"or {DEFAULT_LOGS_DIR})")
    parser.add_argument("--config", type=Path, required=True,
                        help="Path to JSON config file (conda envs, reference paths)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Find and check samples, log what would run, change nothing")
    parser.add_argument("--force", action="store_true",
                        help="Re-run samples even if already marked done")
    parser.add_argument("--max-workers", type=int, default=None,
                        help="Number of samples to process concurrently "
                             "(default: config.json's 'max_parallel_samples', "
                             f"or {DEFAULT_MAX_WORKERS} if that isn't set either)")
    return parser.parse_args(argv)

def main(argv: list[str] | None = None) -> None:
    """Process every sample under the input folder, then run MultiQC.

    Args:
        argv: Command-line arguments; ``None`` uses ``sys.argv``.

    Raises:
        SystemExit: If the config file or input folder is missing, or no
            FASTQ files are found.
    """
    args = parse_args(argv)
    if not args.config.exists():
        raise SystemExit(f"Config file not found: {args.config}")
    with open(args.config) as f:
        config = json.load(f)

    max_workers = max(1, int(args.max_workers or config.get("max_parallel_samples") or DEFAULT_MAX_WORKERS))

    input_dir = args.input.resolve()
    output_dir = args.output.resolve()
    if not input_dir.exists():
        raise SystemExit(f"Input directory does not exist: {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Precedence: --logs-dir > config "logs_dir" > default. One log file per
    # invocation, so a later run on the same day does not overwrite it.
    args.logs_dir = Path(args.logs_dir or config.get("logs_dir", DEFAULT_LOGS_DIR)).resolve()
    args.logs_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.logs_dir / f"pipeline_{datetime.now():%Y%m%d_%H%M%S}.log",
                  with_thread_names=max_workers > 1)

    logger.info("Configuration loaded from %s: %s", args.config, config)
    logger.info("Input: %s  Output: %s  Logs: %s", input_dir, output_dir, args.logs_dir)
    logger.info("Sample-level parallelism: max_workers=%d%s", max_workers, " (DRY RUN)" if args.dry_run else "")

    fastq_files = find_fastq_files(input_dir)
    if not fastq_files:
        raise SystemExit(f"No fastq files found under {input_dir}")

    groups = group_by_sample(fastq_files)
    run_dir = output_dir / f"RESULTS_{datetime.now():%Y%m%d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Discovered %d putative samples; results go to %s", len(groups), run_dir)

    # Threads are enough: the real work runs in `conda run` subprocesses.
    any_completed = False
    with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="sample") as pool:
        futures = {
            pool.submit(process_sample, sample, files, run_dir, input_dir, config, args): sample
            for sample, files in sorted(groups.items())
        }
        for future in as_completed(futures):
            sample = futures[future]
            try:
                any_completed |= future.result()
            except Exception:
                logger.exception("process_sample(%s) raised unexpectedly", sample)

    if args.dry_run or not any_completed:
        return
    multiqc_conda_env = check_env("multiqc_conda_env", config)
    if not multiqc_conda_env:
        logger.warning("multiqc_conda_env not specified or not found. Skipping MultiQC analysis.")
        return
    run_multiqc(multiqc_conda_env, run_dir, run_dir)


if __name__ == "__main__":
    main()