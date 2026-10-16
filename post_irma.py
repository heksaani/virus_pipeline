"""Post-processing steps that run after IRMA.

* ``IRMA_analysis.xlsx``: per-segment N% and read counts for influenza
  (``analyze_and_append_sample``).
* Nextclade: per segment for influenza (``run_nexclade_analysis``) and on the
  whole genome for single-genome species (``run_nextclade_genome_analysis``),
  appended to run-level TSV files.
* MultiQC over the run's FastQC and NanoPlot reports (``run_multiqc``).
* ``sample_summary_report.xlsx``: Kraken2 top hit, ambiguous sites, N% and
  average coverage per sample (``append_sample_row``).

Shared report files are updated under a file lock, so samples running in
parallel can append to them safely.
"""
import fcntl
import logging
import re
from collections.abc import Callable, Generator
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from pre_irma import run_command_capture
from segments import (
    FLU_SEGMENTS,
    consensus_fastas,
    consensus_segments,
    segment_from_filename,
)

logger = logging.getLogger(__name__)

HEADER = [
    "Sample name", "Influenza virus", "Subtype", "HA N%", "NA N%", "PB1 N%",
    "PB2 N%", "PA N%", "MP N%", "NS N%", "NP N%", "Total reads",
    "Matched reads", "Matched reads (%)", "Segment coverage notes",
    "Alternatively matched reads", "Alternatively matched segments",
]

INFL_A_SEGMENT_LENGTHS = {
    "HA": 1701, "NA": 1410, "PB1": 2274, "PB2": 2280,
    "PA": 2151, "MP": 982, "NS": 838, "NP": 1497,
}
INFL_B_SEGMENT_LENGTHS = {
    "HA": 1749, "NA": 1408, "PB1": 2259, "PB2": 2313,
    "PA": 2181, "MP": 1076, "NS": 1027, "NP": 1683,
}

A_SUBTYPE_TO_PREFIX = {
    "H1": "h1n1pdm",
    "H2": "h2n2",
    "H3": "h3n2",
}

# Columns kept from Nextclade's output in the run-level summary TSVs.
NEXTCLADE_SUMMARY_COLUMNS = [
    "seqName", "clade", "subclade", "short-clade", "qc.overallStatus",
    "totalSubstitutions", "totalDeletions", "totalInsertions",
    "totalFrameShifts", "totalMissing", "totalNonACGTNs", "coverage",
    "substitutions", "deletions", "insertions", "frameShifts",
]
NEXTCLADE_SUMMARY_COLUMNS_COVID = [
    "seqName", "clade", "clade_display", "clade_who", "clade_nextstrain",
    "partiallyAliased", "Nextclade_pango", "subclade", "short-clade",
    "qc.overallStatus", "totalSubstitutions", "totalDeletions",
    "totalInsertions", "totalFrameShifts", "totalMissing", "coverage",
    "substitutions", "deletions", "insertions", "frameShifts",
]

# Segment column in nextclade_results.tsv (the lab's Finnish column name).
NEXTCLADE_SEGMENT_COLUMN = "Segmentti"

@contextmanager
def _locked(path: Path) -> Generator[None]:
    """Hold an exclusive lock on ``<path>.lock`` while the block runs.

    The lock is released when the lock file is closed, also if the block
    raises.

    Args:
        path: The shared file to protect.

    Yields:
        Nothing; the lock is held until the block ends.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(Path(f"{path}.lock"), "w") as lockf:
        fcntl.flock(lockf, fcntl.LOCK_EX)
        yield


class Sample:
    """One influenza sample's IRMA results, as a row of ``HEADER``.

    Creating a ``Sample`` never raises because the sample has no output: if
    ``amended_consensus`` is missing or empty, ``valid`` is ``False`` and
    ``result`` keeps its ``None`` placeholders. Check ``valid`` instead of
    wrapping the constructor in try/except.

    Attributes:
        sample_name: Sample name.
        result_path: IRMA results folder that contains ``<sample_name>/``.
        valid: ``True`` if IRMA produced a usable consensus.
        result: One value per ``HEADER`` column.
    """

    def __init__(self, sample_name: str, result_path: Path):
        """Read the sample's IRMA output.

        Args:
            sample_name: Sample name.
            result_path: IRMA results folder that contains ``<sample_name>/``.
        """
        self.sample_name = sample_name
        self.result_path = result_path
        self.valid = self._is_valid()
        self.result: list[Any] = [None] * len(HEADER)
        self.result[0] = self.sample_name
        self.seq_len: dict = {}
        self.coverage_notes = ""

        if self.valid:
            self._get_type()
        if self.valid:
            self._get_segments()
            self._get_read_counts()

    def _is_valid(self) -> bool:
        """Check that IRMA wrote at least one file to ``amended_consensus``.

        Returns:
            ``True`` if ``amended_consensus`` exists and is not empty.
        """
        amended = self.result_path / self.sample_name / "amended_consensus"
        return amended.is_dir() and any(amended.iterdir())

    def _get_type(self) -> None:
        """Fill the type and subtype columns from IRMA's per-gene FASTA names.

        IRMA names them like ``A_HA_H3.fasta`` and ``A_NA_N2.fasta``. For
        influenza B (``B_HA.fasta``) the subtype comes out as ``"HANA"``,
        which ``_normalize_row`` blanks. Marks the sample invalid if there is
        neither an HA nor an NA file.
        """
        ha_type = na_type = infl_type = ""
        for fasta in (self.result_path / self.sample_name).glob("*.fasta"):
            parts = fasta.stem.split("_")
            if "HA" in parts:
                ha_type, infl_type = parts[-1], parts[0]
            elif "NA" in parts:
                na_type, infl_type = parts[-1], parts[0]

        if not ha_type and not na_type:
            logger.warning("Sample %s has neither an HA nor an NA file", self.sample_name)
            self.valid = False
            return

        self.result[1] = infl_type
        self.seq_len = INFL_A_SEGMENT_LENGTHS if infl_type == "A" else INFL_B_SEGMENT_LENGTHS
        self.result[2] = ha_type + na_type

    def _get_segments(self) -> None:
        """Fill the per-segment N% columns and the coverage notes.

        A segment whose consensus is shorter than 95% of its expected length
        is listed in "Segment coverage notes".
        """
        segment_field_map = {
            "HA": 3, "NA": 4, "PB1": 5, "PB2": 6,
            "PA": 7, "MP": 8, "NS": 9, "NP": 10,
        }
        amended = self.result_path / self.sample_name / "amended_consensus"
        for fasta, seg_key in consensus_segments(amended):
            _count, n_percent = get_n_percent_from_fasta(fasta)
            self.result[segment_field_map[seg_key]] = n_percent

            seq_len = _sequence_length(fasta)
            expected = self.seq_len.get(seg_key)
            if expected and seq_len:
                coverage = round(seq_len / expected * 100, 2)
                if coverage < 95:
                    self.coverage_notes += f"{seg_key} coverage {coverage}%, "

        self.result[14] = self.coverage_notes

    def _get_read_counts(self) -> None:
        """Fill the read-count columns from IRMA's ``tables/READ_COUNTS.txt``.

        Counts missing from the file leave their columns empty. If the file
        cannot be read at all, "Total reads" says ``"no results"``.
        """
        reads_txt = self.result_path / self.sample_name / "tables" / "READ_COUNTS.txt"
        try:
            reads_df = pd.read_csv(reads_txt, sep="\t")
        except (OSError, ValueError):  # missing, empty or malformed file
            self.result[11] = "no results"
            return
        if not {"Record", "Reads"} <= set(reads_df.columns):
            self.result[11] = "no results"
            return

        records = reads_df["Record"].astype(str)

        def reads(record: str) -> Any:
            match = reads_df.loc[records == record, "Reads"]
            return match.item() if len(match) == 1 else None

        total, matched = reads("1-initial"), reads("3-match")
        self.result[11] = total
        self.result[12] = matched
        if total and matched is not None:
            self.result[13] = round(matched / total * 100, 2)
        self.result[15] = reads("3-altmatch")
        alt_segments = [r.removeprefix("5-") for r in records if r.startswith("5-")]
        self.result[16] = ", ".join(alt_segments) or None


def _normalize_row(row: list) -> list:
    """Prepare a ``Sample.result`` row for ``IRMA_analysis.xlsx``.

    Empty values become ``"No consensus"`` (empty in the coverage-notes
    column), and the ``"HANA"`` subtype of influenza B becomes empty.

    Args:
        row: One value per ``HEADER`` column.

    Returns:
        The cleaned row.
    """
    coverage_idx = HEADER.index("Segment coverage notes")
    subtype_idx = HEADER.index("Subtype")
    norm = []
    for i, v in enumerate(row):
        if v is None or (isinstance(v, str) and v.strip() == ""):
            norm.append("" if i == coverage_idx else "No consensus")
            continue
        if i == subtype_idx and isinstance(v, str) and v.strip().upper() == "HANA":
            norm.append("")
            continue
        norm.append(v)
    return norm


def _finalize_irma_analysis_df(df: pd.DataFrame) -> pd.DataFrame:
    """Clean up the whole ``IRMA_analysis.xlsx`` table before it is written.

    Blanks the ``"HANA"`` subtype, marks a missing influenza A subtype as
    ``"No consensus"``, and turns the N% columns into numbers (or
    ``"No consensus"``).

    Args:
        df: The combined report table.

    Returns:
        The same table, cleaned.
    """
    if "Subtype" in df.columns:
        df["Subtype"] = df["Subtype"].fillna("").astype(str).str.strip()
        df.loc[df["Subtype"].str.match(r"(?i)^HANA$"), "Subtype"] = ""
        if "Influenza virus" in df.columns:
            mask_empty = df["Subtype"] == ""
            mask_b = df["Influenza virus"] == "B"
            df.loc[mask_empty & ~mask_b, "Subtype"] = "No consensus"

    n_percent_cols = ["HA N%", "NA N%", "PB1 N%", "PB2 N%", "PA N%", "MP N%", "NS N%", "NP N%"]

    def _norm_val(v: Any) -> Any:
        if v is None:
            return "No consensus"
        if isinstance(v, (int, float)):
            return v
        try:
            return float(v)
        except (TypeError, ValueError):
            return "No consensus"

    for c in n_percent_cols:
        if c in df.columns:
            df[c] = df[c].apply(_norm_val)

    return df


def _write_excel(df: pd.DataFrame, excel_path: Path) -> None:
    """Write ``IRMA_analysis.xlsx`` with colour-coded N% and coverage cells.

    Args:
        df: The report table.
        excel_path: Output file.
    """
    excel_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        df.to_excel(writer, sheet_name="N percent", index=False)
        workbook: Any = writer.book
        worksheet: Any = writer.sheets["N percent"]

        zero_fmt = workbook.add_format({"bg_color": "#FFFFFF"})
        worksheet.conditional_format("D2:K1000", {
            "type": "cell", "criteria": "==", "value": 0,
            "format": zero_fmt, "stop_if_true": True,
        })
        worksheet.conditional_format("D2:K1000", {
            "type": "3_color_scale", "min_value": 0, "mid_value": 50, "max_value": 100,
            "max_color": "#003399", "mid_color": "#6699FF", "min_color": "#CCE5FF",
        })
        worksheet.conditional_format("D2:K1000", {
            "type": "text", "criteria": "containing", "value": "No consensus",
            "format": workbook.add_format({"bg_color": "#FF0000"}),
        })
        worksheet.conditional_format("O2:O1000", {
            "type": "text", "criteria": "containing", "value": "coverage",
            "format": workbook.add_format({"bg_color": "#FF0000"}),
        })
        worksheet.set_column("A:A", 25)
        worksheet.set_column("O:O", 35)

def _read_existing_report(path: Path, columns: list[str]) -> pd.DataFrame:
    """Read the current Excel report, or start an empty one.

    An existing report that cannot be read is moved to
    ``<name>.unreadable-<timestamp>`` instead of being overwritten, so its
    rows can be recovered by hand.

    Args:
        path: The report file.
        columns: Columns for a new, empty report.

    Returns:
        The existing rows, or an empty table with ``columns``.

    Raises:
        ImportError: If the Excel reader (openpyxl) is not installed.
    """
    if not path.exists():
        return pd.DataFrame(columns=columns)
    try:
        return pd.read_excel(path)
    except ImportError:
        raise
    except Exception:
        backup = path.with_name(f"{path.name}.unreadable-{datetime.now():%Y%m%d_%H%M%S}")
        path.replace(backup)
        logger.exception("Could not read %s; moved it to %s and started a new report", path, backup)
        return pd.DataFrame(columns=columns)


def _append_row_and_save(row: dict[str, Any] | list,
                         columns: list[str],
                         target_xlsx: Path,
                         dedup_subset: list[str],
                         writer_fn: Callable[[pd.DataFrame, Path], None] | None = None,
                         post_process_fn: Callable[[pd.DataFrame], pd.DataFrame] | None = None) -> None:
    """Add or replace one row in a shared Excel report.

    Safe across samples running in parallel: the file is locked, re-read,
    updated, and replaced through a temporary file.

    Args:
        row: The new row, either a dict keyed by ``columns`` or a list in
            ``columns`` order.
        columns: Column names and order for the report.
        target_xlsx: The Excel file to update.
        dedup_subset: Columns used to drop duplicate rows (keeping the newest)
            and to sort the table, e.g. ``["Sample name"]``.
        writer_fn: How to write the table. Defaults to ``DataFrame.to_excel``;
            ``_write_excel`` adds conditional formatting.
        post_process_fn: Cleanup applied to the whole table before writing.
    """
    with _locked(target_xlsx):
        existing = _read_existing_report(target_xlsx, columns)
        row_df = pd.DataFrame([row], columns=columns)
        combined = pd.concat([existing, row_df], ignore_index=True)
        combined = combined.drop_duplicates(subset=dedup_subset, keep="last")
        combined = combined.sort_values(by=dedup_subset).reset_index(drop=True)

        if post_process_fn:
            combined = post_process_fn(combined)

        tmp_path = target_xlsx.with_suffix(target_xlsx.suffix + ".tmp")
        if writer_fn:
            writer_fn(combined, tmp_path)
        else:
            combined.to_excel(tmp_path, index=False)
        tmp_path.replace(target_xlsx)

    sample_label = row.get(dedup_subset[0]) if isinstance(row, dict) else row[0]
    logger.info("Appended row for %s to %s", sample_label, target_xlsx)


def read_nextclade_tsv(path: Path) -> pd.DataFrame:
    """Read Nextclade's ``nextclade.tsv`` with every value kept as text.

    ``keep_default_na=False`` stops the text "NA" (neuraminidase) from being
    read as a missing value.

    Args:
        path: The ``nextclade.tsv`` file.

    Returns:
        The table, all columns as strings.
    """
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _append_df_to_tsv_locked(df: pd.DataFrame, dest_tsv: Path) -> None:
    """Append rows to a shared, tab-separated Nextclade results file.

    Safe across samples running in parallel. The existing file is read with
    ``keep_default_na=False``, so "NA" segment labels survive each append.

    Args:
        df: Rows to append.
        dest_tsv: The run-level results file.
    """
    with _locked(dest_tsv):
        if dest_tsv.exists():
            existing_df = pd.read_csv(dest_tsv, sep="\t", dtype=str, keep_default_na=False)
            df = pd.concat([existing_df, df], ignore_index=True)
        df.to_csv(dest_tsv, sep="\t", index=False)


def analyze_and_append_sample(result_path: Path, sample_name: str, append_xlsx: Path) -> bool:
    """Add or update a sample's row in the shared ``IRMA_analysis.xlsx``.

    A sample without usable IRMA output still gets a row, filled with
    ``"No consensus"``.

    Args:
        result_path: IRMA results folder that contains ``<sample_name>/``.
        sample_name: Sample name.
        append_xlsx: The ``IRMA_analysis.xlsx`` report.

    Returns:
        ``True`` if IRMA produced a usable consensus for the sample.
    """
    sample = Sample(sample_name, result_path)
    _append_row_and_save(
        row=_normalize_row(sample.result),
        columns=HEADER,
        target_xlsx=append_xlsx,
        dedup_subset=["Sample name"],
        writer_fn=_write_excel,
        post_process_fn=_finalize_irma_analysis_df,
    )
    if not sample.valid:
        logger.warning("Sample %s has no usable amended_consensus output; wrote a 'No consensus' row to %s",
                       sample_name, append_xlsx)
    return sample.valid


def get_dataset_prefix_for_sample(sample_root: Path) -> tuple[str | None, str | None]:
    """Find the influenza type and the Nextclade dataset prefix of a sample.

    Reads IRMA's per-gene FASTA names in the sample folder (not
    ``amended_consensus``), which keep the subtype, e.g. ``A_HA_H1.fasta``.

    Args:
        sample_root: IRMA output folder of the sample.

    Returns:
        ``(infl_type, dataset_prefix)``, e.g. ``("A", "h1n1pdm")`` or
        ``("B", "b")``. The prefix is ``None`` if the H subtype is missing or
        has no dataset; both are ``None`` if the type cannot be determined.
    """
    infl_type = None
    h_token = None

    for fasta in sample_root.glob("*.fasta"):
        parts = fasta.stem.upper().split("_")
        if parts[0] in ("A", "B"):
            infl_type = parts[0]
        if "HA" in parts and len(parts) >= 3 and re.match(r"^H\d+$", parts[-1]):
            h_token = parts[-1]

    if infl_type == "B":
        return infl_type, "b"

    if infl_type == "A":
        if h_token is None:
            logger.warning("No H-subtype token found for sample %s; skipping Nextclade.", sample_root.name)
            return infl_type, None
        prefix = A_SUBTYPE_TO_PREFIX.get(h_token)
        if prefix is None:
            logger.warning("Unmapped H-subtype %s for sample %s; skipping Nextclade.", h_token, sample_root.name)
        return infl_type, prefix

    logger.warning("Could not determine influenza type for sample %s; skipping Nextclade.", sample_root.name)
    return None, None


def run_nexclade_analysis(conda_env: str,
                          nextclade_ref_path: Path,
                          run_folder: Path,
                          input_folder_path: Path | str,
                          output_path: Path,
                          log_path: Path,
                          sample_root: Path,
                          dry_run: bool = False) -> None:
    """Run Nextclade on each influenza segment and collect the results.

    Each segment runs against the dataset
    ``<nextclade_ref_path>/flu_<prefix>_<segment>`` (e.g. ``flu_h3n2_ha``);
    segments without a dataset are skipped. Results are appended to
    ``<run_folder>/nextclade_results.tsv`` with the sample and segment in the
    first two columns.

    Args:
        conda_env: Conda environment with Nextclade.
        nextclade_ref_path: Folder holding the per-segment Nextclade datasets.
        run_folder: Results folder of the run.
        input_folder_path: The sample's ``amended_consensus`` folder.
        output_path: Folder for the per-segment Nextclade output.
        log_path: Log file to append to.
        sample_root: IRMA output folder of the sample; its name is the sample
            name, and its per-gene FASTA names give the subtype.
        dry_run: Log the commands without running them.
    """
    _infl_type, dataset_prefix = get_dataset_prefix_for_sample(sample_root)
    if dataset_prefix is None:
        return

    sample = sample_root.name
    dest_tsv = run_folder / "nextclade_results.tsv"
    for fasta, segment_label in consensus_segments(Path(input_folder_path)):
        segment_ref = nextclade_ref_path / f"flu_{dataset_prefix}_{segment_label.lower()}"
        if not segment_ref.exists():
            logger.warning("Nextclade dataset folder not found: %s, skipping segment %s",
                           segment_ref, segment_label)
            continue

        result = output_path / f"{sample}_{segment_label}_nextclade_results"
        cmd = ["conda", "run", "-n", conda_env, "nextclade", "run",
               "-D", str(segment_ref), "-O", str(result), str(fasta)]
        if dry_run:
            logger.info("DRY RUN Nextclade: %s", " ".join(cmd))
            continue

        logger.info("Running Nextclade for sample %s segment %s (%s)", sample, segment_label, dataset_prefix)
        if not run_command_capture(cmd, log_path):
            logger.error("Nextclade failed for sample %s segment %s", sample, segment_label)
            continue

        nextclade_tsv = result / "nextclade.tsv"
        if not nextclade_tsv.exists():
            logger.warning("Nextclade wrote no nextclade.tsv for sample %s segment %s (%s)",
                           sample, segment_label, result)
            continue

        try:
            df = read_nextclade_tsv(nextclade_tsv)
            df.insert(0, "sample_name", sample)
            df.insert(1, NEXTCLADE_SEGMENT_COLUMN, segment_label)
            cols_to_keep = list(dict.fromkeys(
                ["sample_name", NEXTCLADE_SEGMENT_COLUMN, "seqName"]
                + [c for c in NEXTCLADE_SUMMARY_COLUMNS if c in df.columns]
            ))
            _append_df_to_tsv_locked(df.filter(items=cols_to_keep), dest_tsv)
        except Exception:
            logger.exception("Failed to save Nextclade results for sample %s segment %s",
                             sample, segment_label)
            continue
        logger.info("Nextclade results for sample %s segment %s saved to %s", sample, segment_label, dest_tsv)

def run_nextclade_genome_analysis(conda_env: str,
                                  nextclade_ref_path: Path,
                                  run_folder: Path,
                                  input_folder_path: Path | str,
                                  output_path: Path,
                                  log_path: Path,
                                  dry_run: bool = False,
                                  results_name: str = "nextclade_results_COV.tsv",
                                  summary_columns: list[str] | None = None,
                                  sample: str | None = None) -> None:
    """Run Nextclade on a single-genome sample and collect the result.

    Used for SARS-CoV-2 and any other non-segmented species. The result is
    appended to ``<run_folder>/<results_name>`` with the sample name in the
    first column.

    Args:
        conda_env: Conda environment with Nextclade.
        nextclade_ref_path: The species' Nextclade dataset folder.
        run_folder: Results folder of the run.
        input_folder_path: The sample's ``amended_consensus`` folder.
        output_path: Folder for the per-sample Nextclade output.
        log_path: Log file to append to.
        dry_run: Log the command without running it.
        results_name: File name of the run-level TSV.
        summary_columns: Nextclade columns to keep (all if ``None``).
        sample: Sample name (default: the consensus FASTA's name).
    """
    if not nextclade_ref_path.exists():
        logger.warning("Nextclade dataset folder not found: %s, skipping Nextclade", nextclade_ref_path)
        return
    fasta = next(iter(consensus_fastas(Path(input_folder_path))), None)
    if fasta is None:
        logger.warning("No consensus fasta found in %s", input_folder_path)
        return
    logger.info("Found consensus fasta: %s", fasta)
    sample = sample or fasta.stem

    result = output_path / f"{sample}_nextclade_results"
    cmd = ["conda", "run", "-n", conda_env, "nextclade", "run",
           "-D", str(nextclade_ref_path), "-O", str(result), str(fasta)]
    if dry_run:
        logger.info("DRY RUN Nextclade: %s", " ".join(cmd))
        return

    logger.info("Running Nextclade for sample %s", sample)
    if not run_command_capture(cmd, log_path):
        logger.error("Nextclade failed for sample %s", sample)
        return

    nextclade_tsv = result / "nextclade.tsv"
    if not nextclade_tsv.exists():
        logger.warning("Nextclade wrote no nextclade.tsv for sample %s (%s)", sample, result)
        return

    dest_tsv = run_folder / results_name
    try:
        df = read_nextclade_tsv(nextclade_tsv)
        cols_to_keep = [c for c in (summary_columns or list(df.columns)) if c in df.columns]
        if cols_to_keep:
            df = df.filter(items=cols_to_keep)
        df.insert(0, "sample_name", sample)
        _append_df_to_tsv_locked(df, dest_tsv)
    except Exception:
        logger.exception("Failed to save Nextclade results for sample %s", sample)
        return
    logger.info("Nextclade results for sample %s saved to %s", sample, dest_tsv)


def run_nextclade_covid_analysis(conda_env: str, nextclade_ref_path: Path, run_folder: Path,
                                 input_folder_path: Path | str, output_path: Path, log_path: Path,
                                 dry_run: bool = False) -> None:
    """Run SARS-CoV-2 Nextclade into ``nextclade_results_COV.tsv``.

    Kept for older callers; same as ``run_nextclade_genome_analysis`` with the
    SARS-CoV-2 summary columns.

    Args:
        conda_env: Conda environment with Nextclade.
        nextclade_ref_path: The SARS-CoV-2 Nextclade dataset folder.
        run_folder: Results folder of the run.
        input_folder_path: The sample's ``amended_consensus`` folder.
        output_path: Folder for the per-sample Nextclade output.
        log_path: Log file to append to.
        dry_run: Log the command without running it.
    """
    run_nextclade_genome_analysis(conda_env, nextclade_ref_path, run_folder, input_folder_path,
                                  output_path, log_path, dry_run=dry_run,
                                  results_name="nextclade_results_COV.tsv",
                                  summary_columns=NEXTCLADE_SUMMARY_COLUMNS_COVID)


def run_multiqc(conda_env: str, run_dir: Path, output_dir: Path) -> None:
    """Make one MultiQC report from the run's FastQC and NanoPlot output.

    Looks in ``<run_dir>/<sample>/fastqc_results/`` (Illumina) and
    ``<run_dir>/<sample>/nanoplot/`` (MinION). Skipped if there are none,
    because MultiQC fails without input. ``--force`` overwrites an earlier
    report instead of writing ``multiqc_report_1.html``.

    Args:
        conda_env: Conda environment with MultiQC.
        run_dir: Results folder of the run.
        output_dir: Folder for the report and ``multiqc.log``.
    """
    qc_dirs = [
        str(sample_folder / qc_name)
        for sample_folder in sorted(run_dir.iterdir()) if sample_folder.is_dir()
        for qc_name in ("fastqc_results", "nanoplot")
        if (sample_folder / qc_name).is_dir()
    ]
    if not qc_dirs:
        logger.info("No FastQC/NanoPlot output in %s; skipping MultiQC.", run_dir)
        return

    multiqc_cmd = ["conda", "run", "-n", conda_env, "multiqc", "--force", "-o", str(output_dir), *qc_dirs]
    logger.info("Running MultiQC with command: %s", " ".join(multiqc_cmd))
    if not run_command_capture(multiqc_cmd, output_dir / "multiqc.log"):
        logger.error("MultiQC failed. Check the log at %s", output_dir / "multiqc.log")


def report_header(genome_prefixes: list[str] | tuple[str, ...] = ("Cov",)) -> list[str]:
    """List the columns of ``sample_summary_report.xlsx``.

    Influenza gets three columns per segment; each single-genome species gets
    ``"<prefix> ambiguous count"``, ``"<prefix> N%"`` and
    ``"<prefix> average coverage"``.

    Args:
        genome_prefixes: Column prefixes of the single-genome species, from
            ``species.genome_report_prefixes()``.

    Returns:
        The column names, in order.
    """
    header = ["Sample name", "Kraken species", "Kraken %"]
    for name in (*FLU_SEGMENTS, *genome_prefixes):
        header += [f"{name} ambiguous count", f"{name} N%", f"{name} average coverage"]
    return header


REPORT_HEADER = report_header()  # influenza + SARS-CoV-2 columns, for older callers


def read_report(file_path: str | Path, sep: str = "\t", header: int | None = None) -> pd.DataFrame:
    """Read a tab-separated tool report.

    Args:
        file_path: Report file.
        sep: Column separator.
        header: Header row number, or ``None`` if the file has no header.

    Returns:
        The table, or an empty table if the file is missing or empty.
    """
    path = Path(file_path)
    if not path.exists() or path.stat().st_size == 0:
        logger.warning("Report file %s is missing or empty. Returning empty DataFrame.", path)
        return pd.DataFrame()
    return pd.read_csv(path, sep=sep, header=header)


def get_kraken_information(df: pd.DataFrame) -> tuple[str | None, float | None]:
    """Find the taxon with the most reads assigned directly to it in a Kraken2 report.

    Kraken2 ``--report`` columns (no header): 0 = percentage, 1 = clade reads,
    2 = direct reads, 3 = rank code, 4 = taxid, 5 = name. Every rank is
    considered, including ``unclassified``.

    Args:
        df: The report, as read by ``read_report``.

    Returns:
        ``(name, percentage)``, or ``(None, None)`` for an empty report.
    """
    if df.empty or df.shape[1] < 6:
        return None, None
    top = int(df[2].to_numpy().argmax())
    return str(df[5].to_numpy()[top]).strip(), float(df[0].to_numpy()[top])


def count_low_frequency_consensus(df: pd.DataFrame,
                                  frequency_threshold: float = 0.95) -> tuple[int | None, int | None]:
    """Count consensus positions whose base frequency is below a threshold.

    Args:
        df: IRMA ``*-allAlleles.txt`` table.
        frequency_threshold: Positions with a consensus frequency below this
            count as ambiguous.

    Returns:
        ``(ambiguous_count, total_consensus_positions)``, or ``(None, None)``
        if the table is empty or lacks the needed columns.
    """
    if df.empty or "Frequency" not in df.columns or "Allele_Type" not in df.columns:
        return None, None

    consensus_rows = df[df["Allele_Type"] == "Consensus"]
    total_sites = len(consensus_rows)
    if total_sites == 0:
        return None, None

    low_freq_count = (consensus_rows["Frequency"].astype(float) < frequency_threshold).sum()
    return int(low_freq_count), total_sites

def _segment_files(tables_dir: Path, suffix: str) -> dict[str, Path]:
    """Map each influenza segment to its IRMA table file.

    Args:
        tables_dir: IRMA ``tables`` folder.
        suffix: Table file suffix, e.g. ``"-coverage.txt"`` or
            ``"-allAlleles.txt"``.

    Returns:
        Segment -> file, e.g. ``{"HA": .../A_HA_H3-coverage.txt}``. Files
        whose segment cannot be determined are logged and left out.
    """
    files: dict[str, Path] = {}
    for path in sorted(tables_dir.glob(f"*{suffix}")):
        segment = segment_from_filename(path.name.removesuffix(suffix))
        if segment is None:
            logger.warning("Could not determine segment from filename %s, skipping", path.name)
            continue
        files.setdefault(segment, path)
    return files

def get_genome_ambiguous_sites_from_alleles(alleles_dir: Path) -> int | None:
    """Count ambiguous consensus positions of a single-genome sample.

    Reads the whole-genome ``*-allAlleles.txt`` (e.g.
    ``SARS-CoV-2-allAlleles.txt``); for influenza, which has one file per
    segment, use ``get_flu_ambiguous_sites``.

    Args:
        alleles_dir: IRMA ``tables`` folder.

    Returns:
        The number of ambiguous positions, or ``None`` if there is no usable
        allele table.
    """
    if not alleles_dir.exists():
        logger.warning("Alleles dir not found: %s", alleles_dir)
        return None

    allele_files = sorted(alleles_dir.glob("*-allAlleles.txt"))
    if not allele_files:
        logger.warning("No allAlleles.txt file found in %s", alleles_dir)
        return None
    if len(allele_files) > 1:
        logger.warning("Multiple allAlleles.txt files found in %s, using first: %s",
                       alleles_dir, allele_files[0].name)

    count, _total = count_low_frequency_consensus(read_report(allele_files[0], header=0))
    return count


def get_flu_ambiguous_sites(irma_tables_dir: Path) -> dict[str, tuple[int | None, float | None]]:
    """Count ambiguous consensus positions per influenza segment.

    Args:
        irma_tables_dir: IRMA ``tables`` folder.

    Returns:
        Segment -> ``(ambiguous_count, ambiguous_pct)``, with the percentage
        rounded to two decimals. A segment without consensus positions maps to
        ``(None, None)``; segments without an allele table are left out.
    """
    results: dict[str, tuple[int | None, float | None]] = {}
    if not irma_tables_dir.exists():
        logger.warning("IRMA tables dir not found: %s", irma_tables_dir)
        return results

    for segment, allele_file in _segment_files(irma_tables_dir, "-allAlleles.txt").items():
        count, total = count_low_frequency_consensus(read_report(allele_file, header=0))
        if count is None or not total:
            results[segment] = (None, None)
            continue
        results[segment] = (count, round(count / total * 100, 2))
    return results


def get_genome_ambiguous_sites(snippy_tab_path: Path) -> tuple[int | None, float | None]:
    """Count ambiguous positions in a Snippy-style table (single genome).

    Args:
        snippy_tab_path: The table, with ``Frequency`` and ``Allele_Type``
            columns.

    Returns:
        ``(ambiguous_count, ambiguous_pct)``, or ``(None, None)`` if the file
        is missing or has no consensus positions.
    """
    if not snippy_tab_path.exists():
        logger.warning("Snippy report not found: %s", snippy_tab_path)
        return None, None

    count, total = count_low_frequency_consensus(read_report(snippy_tab_path, header=0))
    if count is None or not total:
        return None, None
    return count, round(count / total * 100, 2)

def _mean_coverage(coverage_file: Path) -> float | None:
    """Return the mean ``Coverage Depth`` of an IRMA coverage table.

    Args:
        coverage_file: IRMA ``*-coverage.txt`` file.

    Returns:
        Mean depth rounded to two decimals, or ``None`` if the table has no
        usable ``Coverage Depth`` column.
    """
    df = read_report(coverage_file, header=0)
    if "Coverage Depth" not in df.columns:
        logger.warning("Coverage column not found in %s", coverage_file)
        return None
    try:
        return round(float(df["Coverage Depth"].astype(float).to_numpy().mean()), 2)
    except (TypeError, ValueError):
        logger.warning("Coverage column in %s is not numeric", coverage_file)
        return None

def get_average_coverage_flu(irma_tables_dir: Path, seq: str) -> float | None:
    """Return the mean coverage depth of one influenza segment.

    Args:
        irma_tables_dir: IRMA ``tables`` folder.
        seq: Segment, e.g. ``"HA"``.

    Returns:
        Mean depth rounded to two decimals, or ``None`` if the segment has no
        usable coverage table.
    """
    coverage_file = _segment_files(irma_tables_dir, "-coverage.txt").get(seq)
    if coverage_file is None:
        logger.warning("Coverage file not found for segment %s in %s", seq, irma_tables_dir)
        return None
    return _mean_coverage(coverage_file)


def get_average_coverage_genome(irma_tables_dir: Path) -> float | None:
    """Return the mean coverage depth of a single-genome sample.

    IRMA names the table after the reference (``<REF>-coverage.txt``), so any
    single-genome species works.

    Args:
        irma_tables_dir: IRMA ``tables`` folder.

    Returns:
        Mean depth rounded to two decimals, or ``None`` if there is no usable
        coverage table.
    """
    coverage_file = next(iter(sorted(irma_tables_dir.glob("*-coverage.txt"))), None)
    if coverage_file is None:
        logger.warning("Coverage file not found in %s", irma_tables_dir)
        return None
    return _mean_coverage(coverage_file)


get_average_coverage_cov = get_average_coverage_genome  # old name, for older callers


def _read_sequence(fasta_path: Path) -> str:
    """Return a FASTA file's sequence lines joined together.

    Args:
        fasta_path: FASTA file.

    Returns:
        The sequence, or ``""`` if the file does not exist.
    """
    if not fasta_path.exists():
        return ""
    with open(fasta_path) as fasta_file:
        return "".join(line.strip() for line in fasta_file if not line.startswith(">"))


def _sequence_length(fasta_path: Path) -> int:
    """Return the sequence length of a FASTA file (0 if it does not exist)."""
    return len(_read_sequence(fasta_path))


def get_n_percent_from_fasta(fasta_path: Path) -> tuple[int | None, float | None]:
    """Count the ``N`` bases in a consensus FASTA.

    Args:
        fasta_path: Consensus FASTA file.

    Returns:
        ``(n_count, n_percent)`` with the percentage rounded to two decimals,
        or ``(None, None)`` if the file is missing or the sequence is empty.
    """
    sequence = _read_sequence(fasta_path)
    if not sequence:
        return None, None
    count = sequence.upper().count("N")
    return count, round(count / len(sequence) * 100, 2)


def get_flu_n_percentages(amended_consensus_dir: Path) -> dict[str, tuple[int | None, float | None]]:
    """Count the ``N`` bases in each influenza segment's consensus.

    Args:
        amended_consensus_dir: IRMA ``amended_consensus`` folder.

    Returns:
        Segment -> ``(n_count, n_percent)``, for segments with a consensus.
    """
    results: dict[str, tuple[int | None, float | None]] = {}
    if not amended_consensus_dir.exists():
        logger.warning("amended_consensus dir not found: %s", amended_consensus_dir)
        return results
    for fasta, segment in consensus_segments(amended_consensus_dir):
        results[segment] = get_n_percent_from_fasta(fasta)
    return results

def get_genome_n_percent(amended_consensus_dir: Path) -> tuple[int | None, float | None]:
    """Count the ``N`` bases in a single-genome sample's consensus.

    Uses the first real (non-padded) consensus FASTA in the folder.

    Args:
        amended_consensus_dir: IRMA ``amended_consensus`` folder.

    Returns:
        ``(n_count, n_percent)``, or ``(None, None)`` if there is no consensus.
    """
    if not amended_consensus_dir.exists():
        logger.warning("amended_consensus dir not found: %s", amended_consensus_dir)
        return None, None
    fasta_file = next(iter(consensus_fastas(amended_consensus_dir)), None)
    if not fasta_file:
        logger.warning("No consensus fasta found in %s", amended_consensus_dir)
        return None, None
    return get_n_percent_from_fasta(fasta_file)


def append_sample_row(sample: str,
                      report_xlsx: Path,
                      segmented: bool = False,
                      genome_prefix: str | None = None,
                      kraken_report_path: Path | None = None,
                      irma_tables_dir: Path | None = None,
                      amended_consensus_dir: Path | None = None,
                      columns: list[str] | None = None) -> bool:
    """Add or update a sample's row in ``sample_summary_report.xlsx``.

    Call once per sample, after all earlier steps have finished. Missing
    inputs leave their columns empty.

    Args:
        sample: Sample name.
        report_xlsx: The summary report file.
        segmented: ``True`` for influenza (per-segment columns).
        genome_prefix: Column prefix for a single-genome species, e.g. ``"Cov"``.
        kraken_report_path: The sample's Kraken2 report.
        irma_tables_dir: The sample's IRMA ``tables`` folder (ambiguous sites
            and coverage).
        amended_consensus_dir: The sample's IRMA ``amended_consensus`` folder
            (N%).
        columns: Report columns (default: ``report_header()`` for this prefix).

    Returns:
        Always ``True``; kept for callers that check it.
    """
    columns = columns or report_header((genome_prefix,) if genome_prefix else ())
    row: dict[str, Any] = dict.fromkeys(columns)
    row["Sample name"] = sample

    if kraken_report_path and kraken_report_path.exists():
        row["Kraken species"], row["Kraken %"] = get_kraken_information(read_report(kraken_report_path))
    else:
        logger.info("No Kraken report for sample %s, leaving species/percentage blank", sample)

    if segmented:
        if irma_tables_dir:
            for seg, (count, _pct) in get_flu_ambiguous_sites(irma_tables_dir).items():
                row[f"{seg} ambiguous count"] = count
                row[f"{seg} average coverage"] = get_average_coverage_flu(irma_tables_dir, seg)
        if amended_consensus_dir:
            for seg, (_count, pct) in get_flu_n_percentages(amended_consensus_dir).items():
                row[f"{seg} N%"] = pct
    elif genome_prefix:
        if irma_tables_dir:
            row[f"{genome_prefix} ambiguous count"] = get_genome_ambiguous_sites_from_alleles(irma_tables_dir)
            row[f"{genome_prefix} average coverage"] = get_average_coverage_genome(irma_tables_dir)
        if amended_consensus_dir:
            row[f"{genome_prefix} N%"] = get_genome_n_percent(amended_consensus_dir)[1]

    _append_row_and_save(row=row, columns=columns, target_xlsx=report_xlsx, dedup_subset=["Sample name"])
    return True