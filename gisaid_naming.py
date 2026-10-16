"""Give influenza samples GISAID names and collect their GISAID submission files.

A THL sample (``YY-NNNNNN`` or ``YYYY-NNNN``) with both an HA and an NA
consensus gets a name such as ``A/Finland/140/2026_H3N2`` or
``B/Finland/141/2026``. The running number comes from the central name list
(``gisaid_naming_file`` in ``config.json``), which is shared by all runs and
maps each THL sample name to its running number and GISAID name.

For each named sample, the run's ``GISAIDfasta/`` folder gets:

* ``<name>/<name>-<SEGMENT>.fasta``: the sample's segments under the new name.
* ``namelist.xlsx``: THL name -> GISAID name.
* ``submit_info.xlsx``: one row per sample for the GISAID batch upload.
* ``fastas.xlsx``: headers and sequences for the upload, one per cell.

All of this happens under one file lock on the central list, so samples
running in parallel never get the same running number.
"""
import fcntl
import logging
import re
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import fastaparser
import pandas as pd
import xlsxwriter

from segments import (
    FLU_SEGMENTS,
    consensus_fastas,
    flu_type_of,
    read_fasta,
    segment_from_filename,
)

logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = ["THL näytenimi", "juokseva numero", "GISAID nimi", "huomoita"]
NAMELIST_COLUMNS = ["THL_name", "GISAID_name"]
SUBMIT_COLUMNS = [
    "Isolate_Id", "", "Isolate_name", "Subtype", "", "", "Location", "", "", "", "Host", "",
    *(f"Seq_Id ({seg})" for seg in FLU_SEGMENTS),
]

# Only THL sample names are GISAID-named: YY-NNNNNN (e.g. 26-057848) or
# YYYY-NNNN. Controls ("Neg-KO-IAV") and MinION barcodes are not.
THL_SAMPLE_NAME = re.compile(r"^(\d{2}|\d{4})-\d+$")


def _has_digit(text: str) -> bool:
    """Return ``True`` if ``text`` contains a digit (``"H3"`` yes, ``"Hx"`` no)."""
    return any(ch.isdigit() for ch in text)


@dataclass
class Segment:
    """One segment's consensus sequence.

    Attributes:
        name: Segment name, e.g. ``"HA"``.
        sequence: Consensus sequence.
    """

    name: str
    sequence: str


class Sample:
    """An influenza sample being GISAID-named.

    Attributes:
        current_name: THL sample name, e.g. ``"26-040001"``.
        running_number: Running number from the central name list.
        sample_path: IRMA output folder of the sample.
        year: Four-digit year from the sample name, or ``"XXXX"``.
        infl_type: ``"A"``, ``"B"`` or ``None``.
        hatype: HA subtype token (``"H3"``), ``""`` if the HA file has none,
            or ``None`` if there is no HA.
        natype: NA subtype token, the same way.
        subtype: Subtype in the name, e.g. ``"H3N2"``, ``"H1N1pdm09"``, ``"B"``.
        new_name: The GISAID name; ``""`` until ``gisaidname()`` succeeds.
        valid: ``False`` if the sample cannot be named.
        error_message: Why the sample is not valid.
        fastas: Consensus FASTA files to submit.
        segments: Segments read from ``fastas``.
        submitting_info: The sample's ``submit_info.xlsx`` row.
    """

    def __init__(self, current_name: str, running_number: str, sample_path: Path):
        """Set up a sample; call ``gisaidname()`` to name it.

        Args:
            current_name: THL sample name, ``YY-NNNNNN`` or ``YYYY-NNNN``.
            running_number: Running number for the GISAID name.
            sample_path: IRMA output folder of the sample.
        """
        self.sample_path = Path(sample_path)
        self.current_name = current_name
        self.running_number = running_number
        self.year = self.get_year()
        self.new_name = ""
        self.segments: list[Segment] = []
        self.fastas: list[Path] = []
        self.infl_type: str | None = None
        self.hatype: str | None = None
        self.natype: str | None = None
        self.subtype: str | None = None
        self.valid = True
        self.submitting_info: list[str] | None = None
        self.error_message = ""

    def get_year(self) -> str:
        """Read the year from the THL sample name.

        ``2026-1234`` gives ``"2026"``; ``26-026817`` gives ``"2026"`` (all
        samples are from the 2000s).

        Returns:
            The four-digit year, or ``"XXXX"`` if the name has no year.
        """
        year_part = self.current_name.split("-")[0]
        if year_part.isdigit() and len(year_part) == 4:
            return year_part
        if year_part.isdigit() and len(year_part) == 2:
            return f"20{year_part}"
        logger.warning("Unrecognized year format in sample name: %s", self.current_name)
        return "XXXX"

    def gisaidname(self) -> None:
        """Build the GISAID name, e.g. ``A/Finland/140/2026_H3N2`` or ``B/Finland/141/2026``.

        Both HA and NA must be present, so no misleading name such as
        ``H1Nx`` is given. Only subtype tokens with a digit are used: if only
        one of HA and NA has one, the name uses that one alone. ``H1N1``
        becomes ``H1N1pdm09``. Sets ``valid`` to ``False`` if the sample
        cannot be named.
        """
        if self.hatype is None or self.natype is None:
            logger.warning("Sample %s is missing HA or NA; not giving it a GISAID name", self.current_name)
            self.valid = False
            self.error_message = "missing HA or NA"
            return

        ha, na = self.hatype, self.natype
        if _has_digit(ha) and _has_digit(na):
            subtype = ha + na
        elif _has_digit(ha):
            subtype = ha
        elif _has_digit(na):
            subtype = na
        else:
            subtype = ""

        base = f"{self.infl_type}/Finland/{self.running_number}/{self.year}"
        if self.infl_type == "A":
            self.subtype = "H1N1pdm09" if subtype == "H1N1" else subtype
            self.new_name = f"{base}_{self.subtype}" if self.subtype else base
        elif self.infl_type == "B":
            self.subtype = "B"
            self.new_name = base
        else:
            logger.warning("Sample %s: influenza type unknown; not giving it a GISAID name", self.current_name)
            self.valid = False
            self.error_message = "unknown influenza type"

    def get_segments(self) -> None:
        """Read each file in ``fastas`` into ``segments``.

        The segment name comes from the file name (see ``segments.py``);
        files whose segment cannot be determined are skipped.
        """
        flu_type = flu_type_of(self.sample_path)
        for fasta in map(Path, self.fastas):
            seg_name = segment_from_filename(fasta.name, flu_type)
            if seg_name is None:
                logger.warning("Cannot tell which segment %s is; not adding it to GISAID files", fasta.name)
                continue
            for header, sequence in read_fasta(fasta):
                logger.info("Adding segment %s (%s)", seg_name, header)
                self.segments.append(Segment(seg_name, sequence))

    def dump_fasta(self, output_path: Path) -> None:
        """Write each segment to ``<output_path>/<name>/<name>-<SEGMENT>.fasta``.

        ``/`` in the GISAID name becomes ``-`` in folder and file names.

        Args:
            output_path: The run's ``GISAIDfasta`` folder.
        """
        name = self.new_name or self.current_name
        fasta_name = name.replace("/", "-")
        fasta_path = Path(output_path) / fasta_name
        fasta_path.mkdir(parents=True, exist_ok=True)
        for segment in self.segments:
            with open(fasta_path / f"{fasta_name}-{segment.name}.fasta", "w") as write_file:
                fastaparser.Writer(write_file).writefasta((f"{name}_{segment.name}", segment.sequence))

    def segments_to_submit(self) -> list[str]:
        """Prepare the sample's GISAID upload entries.

        Also fills ``submitting_info`` with the sample's ``submit_info.xlsx``
        row; its Seq_Id columns are in the order HA, NA, PB1, PB2, PA, MP, NS, NP.

        Returns:
            Alternating FASTA headers and sequences, e.g.
            ``[">A/Finland/1/2026_H3N2_HA", "ACGT...", ...]``.
        """
        seq_ids = dict.fromkeys(FLU_SEGMENTS, "")
        fastas_to_submit = []
        for segment in self.segments:
            seq_id = f"{self.new_name}_{segment.name}"
            fastas_to_submit.extend((f">{seq_id}", segment.sequence))
            if segment.name in seq_ids:
                seq_ids[segment.name] = seq_id

        self.submitting_info = [
            "", "", self.new_name.split("_")[0], self.subtype or "", "", "", "Finland", "", "", "", "human", "",
            *seq_ids.values(),
        ]
        return fastas_to_submit


def is_thl_sample_name(sample_name: str) -> bool:
    """Check whether a sample name is a THL sample name (``YY-NNNNNN`` or ``YYYY-NNNN``).

    Args:
        sample_name: Sample name.

    Returns:
        ``True`` for THL sample names, ``False`` for controls, barcodes, etc.
    """
    return bool(THL_SAMPLE_NAME.match(sample_name))


def detect_type_and_subtype(files: list[Path]) -> tuple[str | None, str | None, str | None]:
    """Find the influenza type and the HA and NA subtypes from FASTA file names.

    IRMA's per-gene files carry both, e.g. ``A_HA_H3.fasta``, ``A_NA_N2.fasta``
    or ``B_HA.fasta``. Name parts are matched exactly (split on ``_``).

    Args:
        files: FASTA files of the sample.

    Returns:
        ``(infl_type, ha, na)``. ``infl_type`` is ``"A"``, ``"B"`` or ``None``.
        ``ha`` is the HA subtype token (``"H3"``), ``""`` if there is an HA file
        without a subtype token, or ``None`` if there is no HA file; ``na`` the
        same for NA.
    """
    infl_type: str | None = None
    subtypes: dict[str, str | None] = {"HA": None, "NA": None}
    for path in sorted(files):
        parts = [p.upper() for p in path.stem.split("_")]
        infl_type = next((p for p in parts if p in ("A", "B")), infl_type)
        for segment, current in subtypes.items():
            if segment not in parts:
                continue
            token = parts[-1]
            if _has_digit(token):
                subtypes[segment] = token
            elif not (current and _has_digit(current)):
                subtypes[segment] = ""
    return infl_type, subtypes["HA"], subtypes["NA"]


@contextmanager
def _file_lock(path: Path) -> Generator[None]:
    """Hold an exclusive lock on ``<path>.lock`` while the block runs.

    Closing the lock file releases the lock, also if the block raises.

    Args:
        path: The shared file to protect.

    Yields:
        Nothing; the lock is held until the block ends.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(Path(f"{path}.lock"), "w") as lockf:
        fcntl.flock(lockf, fcntl.LOCK_EX)
        yield


def _load_names_df(excel_file: Path) -> pd.DataFrame:
    """Read the central name list, or start an empty one.

    Args:
        excel_file: The central name list.

    Returns:
        The list, with all ``REQUIRED_COLUMNS`` present.

    Raises:
        RuntimeError: If the file exists but cannot be read.
    """
    if excel_file.exists():
        try:
            df = pd.read_excel(excel_file)
        except Exception as e:
            raise RuntimeError(f"Error reading central Excel file {excel_file}: {e}") from e
    else:
        df = pd.DataFrame(columns=REQUIRED_COLUMNS)
    for col in REQUIRED_COLUMNS:
        if col not in df.columns:
            df[col] = None
    return df.astype({"GISAID nimi": object, "huomoita": object})


def _read_run_table(path: Path, columns: list[str], match_by_position: bool = False) -> pd.DataFrame:
    """Read one of the run's GISAID tables, or start an empty one.

    A table that cannot be read, or whose columns do not match, is never
    overwritten; naming stops with an error instead, so no earlier rows are
    lost.

    Args:
        path: ``namelist.xlsx`` or ``submit_info.xlsx``.
        columns: Expected columns.
        match_by_position: Rename the columns by position. Needed for
            ``submit_info.xlsx``, whose blank, repeated column names are
            renamed by pandas when the file is read.

    Returns:
        The existing rows, or an empty table with ``columns``.

    Raises:
        RuntimeError: If the file cannot be read or has other columns.
    """
    if not path.exists():
        return pd.DataFrame(columns=columns)
    try:
        df = pd.read_excel(path)
    except Exception as e:
        raise RuntimeError(f"Cannot read {path}: {e}") from e

    if match_by_position and len(df.columns) == len(columns):
        df.columns = columns
    elif match_by_position or not set(columns) <= set(df.columns):
        raise RuntimeError(f"{path} does not have the expected columns; fix or move it and rerun GISAID naming")
    return df


def _read_fasta_list(path: Path) -> list[str]:
    """Read the entries of ``fastas.xlsx`` (one per cell in the first column).

    Args:
        path: ``fastas.xlsx``.

    Returns:
        The entries, or ``[]`` if the file does not exist or is empty.

    Raises:
        RuntimeError: If the file exists but cannot be read.
    """
    if not path.exists():
        return []
    try:
        df = pd.read_excel(path, header=None)
    except Exception as e:
        raise RuntimeError(f"Cannot read {path}: {e}") from e
    return [] if df.empty else df[0].astype(str).tolist()


def _write_run_files(sample: Sample, fastas_to_submit: list[str], output_root: Path) -> None:
    """Add the sample to the run's ``namelist``, ``submit_info`` and ``fastas`` files.

    Args:
        sample: The named sample (``segments_to_submit()`` already called).
        fastas_to_submit: The sample's upload entries.
        output_root: The run's ``GISAIDfasta`` folder.
    """
    namelist_path = output_root / "namelist.xlsx"
    namelist_df = _read_run_table(namelist_path, NAMELIST_COLUMNS)
    new_row = pd.DataFrame([[sample.current_name, sample.new_name]], columns=NAMELIST_COLUMNS)
    namelist_df = pd.concat([namelist_df, new_row], ignore_index=True)
    namelist_df.drop_duplicates(subset=["THL_name"], keep="last").to_excel(namelist_path, index=False)

    submit_path = output_root / "submit_info.xlsx"
    submit_df = _read_run_table(submit_path, SUBMIT_COLUMNS, match_by_position=True)
    submit_df = pd.concat([submit_df, pd.DataFrame([sample.submitting_info], columns=SUBMIT_COLUMNS)],
                          ignore_index=True)
    submit_df["Subtype"] = submit_df["Subtype"].replace("H1N1pdm09", "H1N1")
    submit_df.to_excel(submit_path, index=False)

    fastas_path = output_root / "fastas.xlsx"
    all_fastas = _read_fasta_list(fastas_path) + fastas_to_submit
    workbook = xlsxwriter.Workbook(str(fastas_path))
    worksheet = workbook.add_worksheet()
    for idx, item in enumerate(all_fastas):
        worksheet.write(idx, 0, item)
    workbook.close()


def _running_number(names_df: pd.DataFrame, sample_name: str,
                    running_override: str | None, force_reassign: bool) -> str:
    """Choose the sample's running number.

    Args:
        names_df: The central name list.
        sample_name: THL sample name.
        running_override: Use this number.
        force_reassign: Ignore the number the sample already has.

    Returns:
        ``running_override`` if given; else the sample's existing number
        (unless ``force_reassign``); else the highest number in the list + 1.
    """
    if running_override:
        return str(running_override)
    numbers = pd.Series(pd.to_numeric(names_df["juokseva numero"], errors="coerce"), index=names_df.index)
    existing = numbers.loc[names_df["THL näytenimi"] == sample_name].tolist()
    if not force_reassign and existing and pd.notna(existing[0]):
        return str(int(existing[0]))
    return str(int(numbers.max()) + 1) if numbers.notna().any() else "1"


def _record_name(names_df: pd.DataFrame, sample_name: str, running_number: str,
                 gisaid_name: str) -> pd.DataFrame:
    """Add the sample to the central name list, or update its row.

    Args:
        names_df: The central name list.
        sample_name: THL sample name.
        running_number: The sample's running number.
        gisaid_name: The sample's GISAID name.

    Returns:
        The updated list.
    """
    number = pd.to_numeric(running_number, errors="coerce")
    mask = names_df["THL näytenimi"] == sample_name
    if mask.any():
        idx = names_df.index[mask][0]
        names_df.at[idx, "juokseva numero"] = number
        names_df.at[idx, "GISAID nimi"] = gisaid_name
        names_df.at[idx, "huomoita"] = ""
        return names_df
    row = {"THL näytenimi": sample_name, "juokseva numero": number, "GISAID nimi": gisaid_name, "huomoita": ""}
    return pd.concat([names_df, pd.DataFrame([row])], ignore_index=True)


def process_sample_folder(sample_path: Path,
                          excel_file: Path,
                          run_dir: Path,
                          running_override: str | None = None,
                          force_reassign: bool = False) -> str | None:
    """GISAID-name one sample and write its submission files.

    Non-THL sample names (controls, MinION barcodes) are skipped. The sample
    keeps the running number it already has in the central list unless
    ``force_reassign`` is set.

    Args:
        sample_path: The sample's IRMA output folder (with ``amended_consensus``).
        excel_file: The central THL -> GISAID name list.
        run_dir: Results folder of the run; files go to ``<run_dir>/GISAIDfasta/``.
        running_override: Use this running number.
        force_reassign: Give a new running number even if the sample has one.

    Returns:
        The GISAID name, or ``None`` if the sample was skipped (not a THL name,
        no consensus, or missing HA or NA).

    Raises:
        RuntimeError: If the central list or a run file cannot be read or
            written. Nothing more is changed for the sample then.
    """
    if not is_thl_sample_name(sample_path.name):
        logger.info("Not GISAID-naming %s: not a THL sample name (YY-NNNNNN).", sample_path.name)
        return None

    try:
        with _file_lock(excel_file):
            names_df = _load_names_df(excel_file)
            return _process_sample_locked(sample_path, names_df, excel_file, run_dir,
                                          running_override, force_reassign)
    except OSError as e:
        raise RuntimeError(f"GISAID naming files for {sample_path.name} could not be written: {e}") from e


def _process_sample_locked(sample_path: Path,
                           names_df: pd.DataFrame,
                           excel_file: Path,
                           run_dir: Path,
                           running_override: str | None,
                           force_reassign: bool) -> str | None:
    """Do the work of ``process_sample_folder`` while its lock is held.

    Args:
        sample_path: The sample's IRMA output folder.
        names_df: The central name list, already read.
        excel_file: The central name list's file.
        run_dir: Results folder of the run.
        running_override: Use this running number.
        force_reassign: Give a new running number even if the sample has one.

    Returns:
        The GISAID name, or ``None`` if the sample was skipped.
    """
    sample_name = sample_path.name
    running_number = _running_number(names_df, sample_name, running_override, force_reassign)
    logger.info("Processing sample %s, running number %s", sample_name, running_number)

    if not sample_path.exists():
        logger.error("Sample folder not found: %s", sample_path)
        return None
    amended = sample_path / "amended_consensus"
    if not amended.exists():
        logger.error("amended_consensus folder not found at: %s", amended)
        return None
    fasta_files = consensus_fastas(amended)
    if not fasta_files:
        logger.error("No consensus .fa files found in %s", amended)
        return None
    logger.info("Found %d consensus fasta files", len(fasta_files))

    sample = Sample(sample_name, running_number, sample_path)
    files_to_check = [*sample_path.glob("*.fa*"), *amended.glob("*.fa*")]
    sample.infl_type, sample.hatype, sample.natype = detect_type_and_subtype(files_to_check)
    logger.info("Sample %s: type %s, HA %r, NA %r", sample_name, sample.infl_type, sample.hatype, sample.natype)

    sample.gisaidname()
    if not sample.valid:
        logger.info("Skipping GISAID naming for %s: %s", sample_name, sample.error_message)
        return None
    logger.info("GISAID name: %s", sample.new_name)

    output_root = run_dir / "GISAIDfasta"
    output_root.mkdir(parents=True, exist_ok=True)
    sample.fastas = fasta_files
    sample.get_segments()
    fastas_to_submit = sample.segments_to_submit()
    sample.dump_fasta(output_root)
    _write_run_files(sample, fastas_to_submit, output_root)

    names_df = _record_name(names_df, sample_name, running_number, sample.new_name)
    tmp_path = excel_file.with_suffix(excel_file.suffix + ".tmp")
    names_df.to_excel(tmp_path, index=False)
    tmp_path.replace(excel_file)
    logger.info("GISAID naming complete for %s; files in %s, central list %s updated",
                sample_name, output_root, excel_file)
    return sample.new_name