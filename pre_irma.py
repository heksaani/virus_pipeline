"""Sample discovery and the read-preparation steps that run before IRMA.

Finding and grouping input files:
    find_fastq_files, group_by_sample, get_species_and_seq_type,
    has_r1_r2, unstable_files.

Read preparation per platform (called from workflow.prepare_reads):
    Illumina: FastQC -> fastp -> Kraken2   (run_pre_irma_steps)
    MinION:   merge FASTQs -> NanoPlot -> Kraken2   (run_minion_pre_irma_steps)

Tool wrappers (fastp, FastQC, Kraken2, NanoPlot, Snippy) run each tool in its
own conda environment and append its output to the sample's log file.
"""
import gzip
import logging
import re
import shutil
import subprocess
import time
from functools import lru_cache
from pathlib import Path

logger = logging.getLogger(__name__)

FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")

# Read number at the end of an R1/R2 file name: Illumina/bcl2fastq names
# (X_S1_L001_R1_001.fastq.gz) and simple ones (X_R1.fastq.gz, X_1.fq.gz).
_READ_NUMBER_RE = re.compile(r"[_.\-]R?([12])(?:_\d{3})?\.f(?:ast)?q(?:\.gz)?$", re.IGNORECASE)


def find_fastq_files(input_dir: Path) -> list[Path]:
    """Find all FASTQ files under ``input_dir``, recursively.

    Args:
        input_dir: Folder to search.

    Returns:
        Paths of files ending in ``.fastq``, ``.fq``, ``.fastq.gz`` or ``.fq.gz``.
    """
    files = [p for p in input_dir.rglob("*") if p.is_file() and p.name.endswith(FASTQ_SUFFIXES)]
    logger.info("Found %d FASTQ files under %s", len(files), input_dir)
    return files


@lru_cache(maxsize=1)
def _conda_env_list() -> tuple[tuple[str, str], ...] | None:
    """Return ``(name, path)`` for every conda environment, or ``None`` on error."""
    try:
        result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("Error checking conda environments: %s", e.stderr)
        return None
    except OSError as e:
        logger.error("Could not run conda: %s", e)
        return None

    envs = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # A line looks like "irma_env  *  /opt/conda/envs/irma_env"; the path
        # is always the last column and "*" marks the active environment.
        parts = [p for p in line.split() if p != "*"]
        if parts:
            envs.append((parts[0], parts[-1]))
    return tuple(envs)


def check_conda_env(env_name: str) -> bool:
    """Check whether a conda environment exists, by name or by path.

    A path matches even if it is written differently (relative path,
    symlink), because both sides are resolved before comparing.

    Args:
        env_name: Name or path of the conda environment.

    Returns:
        ``True`` if the environment exists, ``False`` if it does not or the
        environment list could not be read.
    """
    envs = _conda_env_list()
    if envs is None:
        return False
    for name, path in envs:
        if env_name in (name, path):
            return True
        try:
            if Path(env_name).resolve() == Path(path).resolve():
                return True
        except OSError:
            continue
    return False


def check_env(config_key: str, config: dict) -> str | None:
    """Look up a conda environment from the config and check that it exists.

    Logs a warning when the key is missing or the environment does not
    exist, so the caller can simply skip that step.

    Args:
        config_key: Config key holding the environment name, e.g. ``"irma_conda_env"``.
        config: Pipeline configuration.

    Returns:
        The environment name, or ``None`` if the key is missing or the
        environment does not exist.
    """
    env = config.get(config_key)
    if not env:
        logger.warning("%s not specified in config. Skipping step.", config_key)
        return None
    if not check_conda_env(env):
        logger.warning("Conda environment %s (from %s) not found. Skipping step.", env, config_key)
        return None
    return env


def sample_name(filename: str) -> str:
    """Infer the sample name from a FASTQ file name.

    Example: ``26-040001_S12_L001_R1_001.fastq.gz`` -> ``26-040001``.

    Args:
        filename: FASTQ file name or path.

    Returns:
        The part before ``_S<digits>`` or ``_R1``/``_R2``, or the file name
        up to its first dot if neither is found.
    """
    bn = Path(filename).name
    m = re.match(r"(?P<s>.+?)_S\d+", bn) or re.match(r"(?P<s>.+?)_R[12]", bn)
    if m:
        return m.group("s")
    return bn.split(".")[0]


def is_file_stable(path: Path, wait_seconds: float = 3.0) -> bool:
    """Check whether a file has finished being written.

    A file can show up in a directory listing while it is still being copied.
    Comparing its size before and after a short wait catches that.

    Args:
        path: File to check.
        wait_seconds: Time between the two size checks.

    Returns:
        ``True`` if the file exists, is not empty and did not grow during
        the wait.
    """
    return not unstable_files([path], wait_seconds)


def unstable_files(paths: list[Path], wait_seconds: float = 3.0) -> list[Path]:
    """Find files that are missing, empty or still growing.

    All files are checked with a single wait, so a MinION barcode with
    hundreds of files costs one ``wait_seconds``, not one per file.

    Args:
        paths: Files to check.
        wait_seconds: Time between the two size checks (0 in tests).

    Returns:
        The files from ``paths`` that are missing, empty or changed size
        during the wait; an empty list if all are stable.
    """
    def size(p: Path) -> int | None:
        try:
            return p.stat().st_size
        except FileNotFoundError:
            return None

    before = {p: size(p) for p in paths}
    if wait_seconds > 0:
        time.sleep(wait_seconds)
    return [p for p in paths if not before[p] or size(p) != before[p]]

def group_by_sample(files: list[Path]) -> dict[str, list[Path]]:
    """Group FASTQ files by sample name.

    A file inside a ``barcode*`` folder (MinION) belongs to that barcode.
    Files in MinKNOW's ``unclassified`` folder are grouped as one
    ``"unclassified"`` sample, which ``run_pipeline`` skips. Any other file
    is grouped by ``sample_name()`` of its file name.

    Args:
        files: FASTQ files.

    Returns:
        Sample name -> files of that sample.
    """
    groups: dict[str, list[Path]] = {}
    for p in files:
        barcode_dir = next(
            (parent.name for parent in p.parents
             if parent.name.lower().startswith("barcode") or parent.name.lower() == "unclassified"),
            None,
        )
        groups.setdefault(barcode_dir or sample_name(p.name), []).append(p)
    return groups


def get_species_and_seq_type(files_sorted: list[str], input_dir: Path) -> tuple[str | None, str | None]:
    """Find the species code and sequencer type from a sample's file path.

    Expected layout::

        <input_dir>/<species folder>/<sequencer folder>/.../file.fastq.gz

    The species folder is matched against ``species.SPECIES`` (folder names
    and aliases); the sequencer folder against ``species.SEQUENCER_PLATFORMS``.

    Args:
        files_sorted: The sample's files; the first one is used.
        input_dir: Pipeline input folder.

    Returns:
        ``(species_code, seq_type)``, e.g. ``("FLU", "nextseq")``. Either part
        is ``None`` if it could not be determined.
    """
    if not files_sorted:
        logger.warning("No files provided to determine species/seq_type")
        return None, None

    file_path = Path(files_sorted[0])
    try:
        rel_parts = file_path.resolve().relative_to(input_dir.resolve()).parts
    except ValueError:
        logger.warning("Path %s is not under input_dir %s", file_path, input_dir)
        return None, None

    if len(rel_parts) < 2:
        logger.warning("Path %s does not have expected <species>/<seq_type>/... structure", file_path)
        return None, None

    species_folder, seq_type_folder = rel_parts[0].lower(), rel_parts[1].lower()

    # Imported here: species.py imports the workflows, which import this module.
    from species import SEQUENCER_PLATFORMS, species_by_folder

    species_def = species_by_folder(species_folder)
    species = species_def.code if species_def else None
    if species is None:
        logger.warning("Unknown species folder '%s' in path %s", species_folder, file_path)

    seq_type = seq_type_folder if seq_type_folder in SEQUENCER_PLATFORMS else None
    if seq_type is None:
        logger.warning("Unknown sequencer folder '%s' in path %s", seq_type_folder, file_path)

    logger.info("Determined species: %s, sequencer type: %s", species, seq_type)
    return species, seq_type

def read_number(filename: str) -> int | None:
    """Read the R1/R2 read number from a FASTQ file name.

    Only the end of the name is checked, so a sample name containing
    ``r1`` or ``_2_`` is not mistaken for a read number.

    Args:
        filename: FASTQ file name.

    Returns:
        1 or 2, or ``None`` if the name is not an R1/R2 file.
    """
    m = _READ_NUMBER_RE.search(filename)
    return int(m.group(1)) if m else None


def has_r1_r2(file_paths: list[str]) -> tuple[bool, bool]:
    """Check which reads of a pair are present.

    Args:
        file_paths: FASTQ file paths.

    Returns:
        ``(has_r1, has_r2)``.
    """
    reads = {read_number(Path(fp).name) for fp in file_paths}
    return 1 in reads, 2 in reads

def combine_fastqs(files: list[Path | str], dest: Path) -> Path:
    """Concatenate FASTQ files (gzipped and/or plain) into one gzipped file.

    Gzipped inputs are copied as-is (concatenated gzip members are valid
    gzip); plain inputs are compressed on the way. Data is streamed in
    chunks, and the result is written to ``<dest>.part`` first, so an
    interrupted run never leaves a truncated file that looks complete.

    Args:
        files: FASTQ files to combine, in any order (they are sorted).
        dest: Output ``.fastq.gz`` path.

    Returns:
        ``dest``.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    with open(tmp, "wb") as out:
        for f in sorted(map(Path, files)):
            if f.name.endswith(".gz"):
                with open(f, "rb") as src:
                    shutil.copyfileobj(src, out, 1024 * 1024)
            else:
                with open(f, "rb") as src, gzip.GzipFile(fileobj=out, mode="wb") as gz:
                    shutil.copyfileobj(src, gz, 1024 * 1024)
    tmp.replace(dest)
    return dest

def translate_paths_for_container(files_sorted: list[str],
                                  input_dir: Path, run_dir: Path) -> list[str]:
    """Map host paths to the paths a Docker container sees.

    ``input_dir`` is mounted at ``/data/input`` and ``run_dir`` at
    ``/data/output``.

    Args:
        files_sorted: Host paths to translate.
        input_dir: Host folder mounted as ``/data/input``.
        run_dir: Host folder mounted as ``/data/output``.

    Returns:
        Container paths, in the same order. Paths outside both folders are
        returned unchanged.
    """
    mounts = ((input_dir.resolve(), "/data/input"), (run_dir.resolve(), "/data/output"))
    mapped = []
    for p in files_sorted:
        resolved = Path(p).resolve()
        target = next((f"{mount}/{resolved.relative_to(host)}"
                       for host, mount in mounts if resolved.is_relative_to(host)), str(Path(p)))
        mapped.append(target)
    return mapped


def run_command(cmd: list[str], cwd: Path | None = None) -> None:
    """Run a command and log it.

    Args:
        cmd: Command and arguments.
        cwd: Working folder (optional).

    Raises:
         If the command exits with a non-zero status.
    """
    logger.info("Running: %s", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True, cwd=cwd)


def run_command_capture(cmd: list[str], log_path: Path, cwd: Path | None = None) -> bool:
    """Run a command and append its stdout and stderr to a log file.

    Args:
        cmd: Command and arguments.
        log_path: Log file to append to.
        cwd: Working folder (optional).

    Returns:
        ``True`` if the command exited with 0; ``False`` if it failed or
        could not be started (e.g. binary not found).
    """
    logger.info("Running (captured): %s -> %s", " ".join(map(str, cmd)), log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write(f"\n--- COMMAND: {' '.join(map(str, cmd))} ---\n")
        fh.flush()
        try:
            subprocess.run(cmd, check=True, cwd=cwd, stdout=fh, stderr=fh)
            return True
        except subprocess.CalledProcessError as e:
            fh.write(f"\n--- COMMAND FAILED (exit code {e.returncode}) ---\n")
            logger.error("Command failed (exit code %s): %s", e.returncode, " ".join(map(str, cmd)))
            return False
        except OSError as e:
            fh.write(f"\n--- COMMAND FAILED TO START: {e} ---\n")
            logger.error("Command failed to start: %s (%s)", " ".join(map(str, cmd)), e)
            return False


def run_fastp(
    conda_env: str,
    sample: str,
    files: list[str],
    out_dir: Path,
    log_path: Path,
    dry_run: bool = False,
    threads: int = 4,
) -> list[str]:
    """Trim and filter reads with fastp.

    Args:
        conda_env: Conda environment with fastp.
        sample: Sample name, used in output file names.
        files: One file (single-end) or two (R1, R2). Extra files are ignored
            with a warning.
        out_dir: Output folder for trimmed reads and the HTML/JSON reports.
        log_path: Log file to append to.
        dry_run: Log the command without running it.
        threads: Number of threads.

    Returns:
        Trimmed FASTQ paths (one or two), or ``[]`` if there was no input or
        fastp failed.
    """
    if not files:
        logger.warning("No FASTQ files provided for fastp, skipping")
        return []
    if len(files) > 2:
        logger.warning("More than 2 FASTQ files provided for fastp, only using the first two")
        files = files[:2]

    base_cmd = ["conda", "run", "-n", conda_env, "fastp", "--thread", str(threads),
                "--html", str(out_dir / f"{sample}_fastp_report.html"),
                "--json", str(out_dir / f"{sample}_fastp_report.json")]

    if len(files) == 1:
        logger.info("Single-end FASTQ file provided for fastp: %s", files[0])
        output = str(out_dir / f"{sample}_fastp.fastq.gz")
        cmd = base_cmd + ["-i", files[0], "-o", output]
        results = [output]
    else:
        logger.info("Paired-end FASTQ files provided")
        output_r1 = str(out_dir / f"{sample}_fastp_R1.fastq.gz")
        output_r2 = str(out_dir / f"{sample}_fastp_R2.fastq.gz")
        cmd = base_cmd + ["-i", files[0], "-o", output_r1, "-I", files[1], "-O", output_r2]
        results = [output_r1, output_r2]
    if dry_run:
        logger.info("DRY RUN fastp: %s", " ".join(cmd))
        return results
    logger.info("Running fastp for sample %s -> %s", sample, out_dir)

    if not run_command_capture(cmd, log_path):
        logger.error("fastp failed for sample %s", sample)
        return []
    return results


def run_fastqc(
    conda_env: str | None,
    files: list[str],
    out_dir: Path,
    log_path: Path,
    dry_run: bool = False,
    threads: int = 2,
) -> None:
    """Write FastQC reports for FASTQ files.

    Args:
        conda_env: Conda environment with FastQC, or ``None`` to skip this step.
        files: FASTQ files.
        out_dir: Output folder for the reports.
        log_path: Log file to append to.
        dry_run: Log the command without running it.
        threads: Number of threads.
    """
    if not conda_env:
        logger.warning("No conda env provided for FastQC; skipping.")
        return
    cmd = [
        "conda", "run", "-n", conda_env,
        "fastqc",
        "--outdir", str(out_dir),
        "--threads", str(threads),
        "--quiet",
        *files,
    ]
    if dry_run:
        logger.info("DRY RUN FastQC: %s", " ".join(cmd))
        return
    logger.info("Running FastQC for %d file(s): %s -> %s", len(files), ", ".join(files), out_dir)
    if not run_command_capture(cmd, log_path):
        logger.error("FastQC failed for files: %s", files)

def run_snippy(
    conda_env: str,
    r1: str,
    r2: str,
    ref: Path,
    out_dir: Path,
    log_path: Path,
    dry_run: bool = False,
    threads: int = 4,
) -> None:
    """Call variants against a reference with Snippy (paired-end reads).

    ``run_pipeline`` reads ``snps.tab`` from ``run_dir/<sample>/snippy/``,
    so callers pass that folder as ``out_dir``.

    Args:
        conda_env: Conda environment with Snippy.
        r1: R1 FASTQ path.
        r2: R2 FASTQ path.
        ref: Reference genome (``.gbk`` or ``.fa``).
        out_dir: Folder Snippy writes into (overwritten with ``--force``).
        log_path: Log file to append to.
        dry_run: Log the command without running it.
        threads: Number of threads.
    """
    cmd = [
        "conda", "run", "-n", conda_env,
        "snippy",
        "--cpus", str(threads),
        "--outdir", str(out_dir),
        "--ref", str(ref),
        "--R1", r1,
        "--R2", r2,
        "--force",
    ]
    if dry_run:
        logger.info("DRY RUN Snippy: %s", " ".join(cmd))
        return
    logger.info("Running Snippy -> %s", out_dir)
    if not run_command_capture(cmd, log_path):
        logger.error("Snippy failed (output folder %s)", out_dir)

def run_kraken(
    conda_env: str,
    reads: list[str],
    db: Path | str,
    report_path: Path,
    output_path: Path,
    log_path: Path,
    dry_run: bool = False,
    threads: int = 4,
) -> None:
    """Classify reads taxonomically with Kraken2.

    A read pair is run with ``--paired``, so mates are classified together.

    Args:
        conda_env: Conda environment with Kraken2.
        reads: ``[reads.fastq.gz]`` (MinION) or ``[R1, R2]`` (Illumina).
        db: Kraken2 database folder.
        report_path: Output path for the Kraken2 report.
        output_path: Output path for the per-read classifications.
        log_path: Log file to append to.
        dry_run: Log the command without running it.
        threads: Number of threads.
    """
    cmd = ["conda", "run", "-n", conda_env, "kraken2",
           "--db", str(db),
           "--threads", str(threads),
           "--report", str(report_path),
           "--output", str(output_path)]
    if len(reads) == 2:
        cmd.append("--paired")
    cmd += [str(r) for r in reads]
    if dry_run:
        logger.info("DRY RUN Kraken: %s", " ".join(cmd))
        return
    logger.info("Running Kraken for %s -> %s", ", ".join(map(str, reads)), report_path)
    if not run_command_capture(cmd, log_path):
        logger.error("Kraken failed for %s", ", ".join(map(str, reads)))

def run_kraken_if_configured(sample: str, reads: list[str], sample_dir: Path,
                             log_path: Path, config: dict) -> None:
    """Run Kraken2 into ``<sample_dir>/kraken/`` if its env and database are set.

    Kraken2 is a QC report, not an input to IRMA, so a missing environment
    or database only skips this step.

    Args:
        sample: Sample name, used in output file names.
        reads: Reads to classify (see ``run_kraken``).
        sample_dir: Results folder for this sample.
        log_path: Log file to append to.
        config: Pipeline configuration (``kraken_conda_env``, ``kraken-db``,
            ``kraken_threads``).
    """
    kraken2_conda_env = check_env("kraken_conda_env", config)
    kraken_db_path = config.get("kraken-db")
    if not (kraken2_conda_env and kraken_db_path):
        logger.warning("kraken_conda_env or kraken-db missing; skipping Kraken2 for %s.", sample)
        return
    kraken_dir = sample_dir / "kraken"
    kraken_dir.mkdir(parents=True, exist_ok=True)
    run_kraken(kraken2_conda_env, reads, kraken_db_path,
               kraken_dir / f"{sample}_kraken_report.txt",
               kraken_dir / f"{sample}_kraken_output.txt",
               log_path=log_path, threads=config.get("kraken_threads", 4))

def run_pre_irma_steps(sample: str, files_sorted: list, sample_dir: Path,
                       log_path: Path, config: dict) -> list[str] | None:
    """Prepare one paired-end Illumina sample: FastQC -> fastp -> Kraken2.

    FastQC and Kraken2 are optional; a missing environment only logs a
    warning. fastp is required because IRMA runs on its output.

    Args:
        sample: Sample name.
        files_sorted: The sample's raw FASTQ files (R1 and R2).
        sample_dir: Results folder for this sample.
        log_path: Log file to append to.
        config: Pipeline configuration.

    Returns:
        ``[trimmed_R1, trimmed_R2]``, or ``None`` if R1 or R2 is missing or
        fastp could not run.
    """
    r1, r2 = has_r1_r2(files_sorted)
    if not r1 or not r2:
        logger.warning("Sample %s does not have both R1 and R2 files. Skipping.", sample)
        return None

    fastqc_out = sample_dir / "fastqc_results"
    fastqc_out.mkdir(parents=True, exist_ok=True)
    run_fastqc(check_env("qc_conda_env", config), files_sorted, fastqc_out, log_path,
               threads=config.get("fastqc_threads", 4))

    fastp_conda_env = check_env("fastp_conda_env", config)
    if not fastp_conda_env:
        return None
    fastp_results = run_fastp(fastp_conda_env, sample, files_sorted, sample_dir, log_path,
                              threads=config.get("fastp_threads", 4))
    if len(fastp_results) != 2:
        logger.error("fastp did not produce a read pair for sample %s (got %s). Skipping further steps.",
                     sample, fastp_results)
        return None

    run_kraken_if_configured(sample, fastp_results, sample_dir, log_path, config)
    return fastp_results


def run_nanoplot(conda_env: str, reads: Path, out_dir: Path, log_path: Path, threads: int = 4) -> None:
    """Write a NanoPlot read-length and quality report for MinION reads.

    MultiQC picks the report up. (FastQC and fastp are Illumina tools.)

    Args:
        conda_env: Conda environment with NanoPlot.
        reads: Merged reads file.
        out_dir: Output folder for the report.
        log_path: Log file to append to.
        threads: Number of threads.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", conda_env, "NanoPlot", "--fastq", str(reads),
           "--outdir", str(out_dir), "--threads", str(threads)]
    if not run_command_capture(cmd, log_path):
        logger.error("NanoPlot failed for %s", reads)

def run_minion_pre_irma_steps(sample: str, files_sorted: list, sample_dir: Path,
                              log_path: Path, config: dict) -> Path | None:
    """Prepare one MinION barcode: merge reads -> NanoPlot -> Kraken2.

    All of the barcode's FASTQs are used, from both ``fastq_pass`` and
    ``fastq_fail`` (lab decision 2026-10-06). Reads are not trimmed, because
    IRMA's MinION configs are built for raw ONT reads. NanoPlot and Kraken2
    are optional.

    Args:
        sample: Sample (barcode) name.
        files_sorted: The barcode's FASTQ files.
        sample_dir: Results folder for this sample.
        log_path: Log file to append to.
        config: Pipeline configuration.

    Returns:
        Path to ``<sample_dir>/<sample>.fastq.gz``, or ``None`` if there were
        no files.
    """
    if not files_sorted:
        logger.warning("Sample %s has no FASTQ files.", sample)
        return None
    reads = combine_fastqs(files_sorted, sample_dir / f"{sample}.fastq.gz")
    logger.info("Merged %d FASTQ files for %s into %s", len(files_sorted), sample, reads)

    nanoplot_env = check_env("nanoplot_conda_env", config)
    if nanoplot_env:
        run_nanoplot(nanoplot_env, reads, sample_dir / "nanoplot", log_path,
                     threads=config.get("nanoplot_threads", 4))

    run_kraken_if_configured(sample, [str(reads)], sample_dir, log_path, config)
    return reads