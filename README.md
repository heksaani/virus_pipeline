# IRMA Pipeline Runner

A pipeline for influenza and SARS-CoV-2 sequencing data from Illumina
(NextSeq / MiSeq) and Oxford Nanopore MinION runs. It runs read QC,
taxonomic classification, assembly with IRMA, clade assignment and GISAID
naming, and writes Excel/HTML summary reports.

| Species | NextSeq | MiSeq | MinION |
|---|---|---|---|
| Influenza (`influenza/`) | ✅ | ✅ | ✅ (no GISAID naming) |
| SARS-CoV-2 (`covid/`) | ✅ | ✅ | ✅ (no Snippy) |

## Pipeline steps

1. **Checks** – files no longer growing; for Illumina, both R1 and R2 present.
2. **Pre-processing**
   - Illumina: `FastQC` (raw reads), `fastp` (trimming), `Kraken2` (trimmed pairs).
   - MinION: all of the barcode's reads (`fastq_pass` **and** `fastq_fail`,
     gzipped or plain) merged into `<sample>.fastq.gz`, `NanoPlot`, `Kraken2`.
     No trimming or filtering.
3. **Assembly** – `IRMA` with the module set for the species and platform in
   `species.py`: `FLU` / `CoV` for Illumina, `FLU-minion` /
   `CoV-minion-long-reads` for MinION; plus `Snippy` for Illumina SARS-CoV-2.
4. **Post-processing**
   - Influenza: consensus segments copied to `segments/<SEGMENT>/`,
     `IRMA_analysis.xlsx` (N% per segment, read counts), GISAID naming
     (Illumina only), `Nextclade` per segment (A/H1, A/H2, A/H3 and B datasets).
   - SARS-CoV-2: `Nextclade` on the consensus genome.
   - All: `sample_summary_report.xlsx` (Kraken top hit, ambiguous sites,
     N%, mean coverage) and one `MultiQC` report per run.
5. **Finish** – raw FASTQs moved to `archive_dir`, sample marked done.

## Requirements

- Python 3.11+ (developed on 3.13) and [uv](https://docs.astral.sh/uv/)
- conda, with each tool in an environment named in `config.json`
  (fastp, FastQC, Kraken2, IRMA, Snippy, Nextclade, MultiQC; NanoPlot for MinION)
- Reference data: a Kraken2 database, the SARS-CoV-2 reference (`.gbk`) for
  Snippy, and Nextclade datasets (see [Nextclade datasets](#nextclade-datasets))

## Installation

```bash
git clone <repository_url>
cd irma_pipe
uv lock      # regenerate uv.lock after the 2026-10-06 pyproject.toml change (adds openpyxl)
uv sync      # installs fastaparser, pandas, openpyxl, xlsxwriter, watchdog (+ pytest, ruff for dev)

for env_file in conda_envs/*.yml; do conda env create -f "$env_file"; done
```

> `openpyxl` is required: pandas needs it to read the existing `.xlsx`
> reports and the central GISAID name list.

## Configuration (`config.json`)

| Key | Used for | Required |
|---|---|---|
| `qc_conda_env` | FastQC | no (step skipped) |
| `fastp_conda_env` | fastp | **yes** (Illumina) |
| `nanoplot_conda_env` | NanoPlot read report (MinION) | no (step skipped) |
| `kraken_conda_env`, `kraken-db` | Kraken2 | no (step skipped) |
| `irma_conda_env` | IRMA | **yes** |
| `irma_external_config` | IRMA `--external-config` file for all species | no – defaults to `custom-config.sh` next to the code |
| `irma_external_config_<code>` | Per-species IRMA config, e.g. `irma_external_config_flu` | no – defaults to `custom-config-<code>.sh` next to the code if it exists |
| `irma_module_<code>_<platform>` | Override a species' IRMA module, e.g. `irma_module_cov_minion` | no – defaults in `species.py` |
| `snippy_conda_env`, `sars-ref` | Snippy (CoV) | no |
| `nextclade_conda_env`, `nextclade_ref` | Nextclade, influenza | no |
| `nextclade_ref_covid` | Nextclade, SARS-CoV-2 | no |
| `multiqc_conda_env` | MultiQC | no |
| `gisaid_naming_file` | Central THL → GISAID name list (influenza) | no (naming skipped) |
| `logs_dir` | Per-run and per-sample logs, `.IRMA_DONE` markers | default `/mnt/volume/logs` |
| `archive_dir` | Where raw FASTQs go after a sample completes | default `/mnt/volume/storage` |
| `max_parallel_samples` | Samples processed at once | default 2 |
| `fastqc_threads`, `fastp_threads`, `kraken_threads`, `snippy_threads`, `nanoplot_threads` | Threads per tool | default 4 |
| `file_stable_wait_seconds` | How long to watch a FASTQ's size before trusting it | default 3 |

Required steps stop the sample (raw files stay in the input folder for a
retry). Optional steps log a warning and the sample continues.

IRMA reference: <https://github.com/CDCgov/irma> (modules and their configs
are in `IRMA_RES/modules/<MODULE>/`).

`custom-config.sh` sets IRMA options (e.g. `MIN_CONS_SUPPORT=20`) and IRMA's
`TMP` folder. Check that `TMP` points at a folder that exists on your server.
Influenza uses `custom-config-flu.sh` instead: the same settings plus the
lab's `SEG_NUMBERS` (see below). Keep the shared lines of the two files in
step. IRMA applies the external config after the module's own `init.sh` and
config, so these values win.

### Nextclade datasets

`nextclade_ref` must contain one folder per subtype and segment, named
`flu_<prefix>_<segment>`, e.g. `flu_h3n2_ha`, `flu_h1n1pdm_na`, `flu_b_ha`.
Segments without a matching folder are skipped with a warning.

## Expected input layout

```
<input_dir>/<influenza|covid>/<nextseq|miseq|minion>/<RUN_FOLDER>/.../<sample>_S<n>_L001_R1_001.fastq.gz
```

Folder names are case-insensitive. **Copy each sequencing run in as its own
folder** – the results are named after it. A MiSeq run can be copied in as-is:

```
input/influenza/miseq/261006_M01234_0042_000000000-ABCDE/Data/Intensities/BaseCalls/
    26-040001_S1_L001_R1_001.fastq.gz
    26-040001_S1_L001_R2_001.fastq.gz
    Undetermined_S0_L001_R1_001.fastq.gz      <- skipped, left in place
```

Results for this run go to `<results>/261006_M01234_0042_000000000-ABCDE_RESULTS/`,
however many pipeline runs it takes and whatever the time of day.
FASTQs placed directly in `<species>/<seq_type>/` (no run folder) still work:
their results are named after the date of the sample's oldest FASTQ
(`20261006_RESULTS`), and a warning is logged.

The sample name is the part before `_S<number>` (here `26-040001`); it is
also the THL sample name used for GISAID naming, so it must be in
`YY-NNNNNN` or `YYYY-NNNN` form. `Undetermined` reads are never analysed.

**MinION:** copy the finished run folder (the one MinKNOW writes, with
`fastq_pass/` and `fastq_fail/`) under `<species>/minion/`:

```
input/influenza/minion/20261006_1203_MN41234_FAX12345_abcd1234/
    fastq_pass/barcode01/*.fastq.gz      <- sample "barcode01"
    fastq_fail/barcode01/*.fastq.gz      <- same sample, also used
    fastq_pass/unclassified/             <- skipped, left in place
```

Each `barcodeNN` folder is one sample, named `barcodeNN` throughout the
results. Because these are not THL sample names, MinION samples get no
GISAID name. The pipeline does not wait for MinKNOW to finish, so only copy
runs that have completed.

## Usage

### One-off run

```bash
uv run python run_pipeline.py -i /path/to/input -o /path/to/results --config config.json
```

| Option | Meaning |
|---|---|
| `--config` | **required** JSON config |
| `--logs-dir` | overrides `logs_dir` from the config |
| `--dry-run` | find and check samples, log what would run; runs no tools, moves no files |
| `--force` | re-run samples that already have a `.IRMA_DONE` marker |
| `--max-workers N` | samples processed at once (overrides `max_parallel_samples`) |

Each sequencing run gets `<results>/<RUN_FOLDER>_RESULTS/`, reused by every
later invocation that sees files from that run (late R2, retries, watcher
re-runs), including after midnight. Each invocation writes its own log,
`<logs_dir>/pipeline_<YYYYMMDD_HHMMSS>.log`. Each sample has a
`pipeline.log` with full tool output next to its `.IRMA_DONE` marker in
`<logs_dir>/<RUN_FOLDER>_RESULTS/<sample>/`. Raw FASTQs are archived to
`<archive_dir>/<RUN_FOLDER>_RESULTS/<sample>/`.

### What happens to the raw files

| Situation | Raw FASTQs | `.IRMA_DONE` |
|---|---|---|
| Sample completed (including "IRMA found no consensus") | moved to `archive_dir/<run>/<sample>/` (MinION: `…/<sample>/fastq_pass/` and `fastq_fail/`) | written |
| Only R1 or R2 has arrived, or a file is still growing | left in place | – |
| A required step failed (fastp, IRMA, missing conda env) | left in place | – |
| Unexpected error | left in place, traceback in the log | – |
| Unknown species/sequencer folder, `Undetermined` or MinION `unclassified` | left in place | – |
| `--dry-run` | left in place | – |

Anything left in place is picked up by the next run, so fix the cause and
run again (or let the watcher do it).

### Continuous / watch mode

```bash
uv run python watch_and_run.py -i /path/to/input -o /path/to/results --config config.json --debounce 60
```

Watches `input_dir` recursively. After `--debounce` seconds with no new
FASTQ files it runs `run_pipeline.py` (next to the watcher, override with
`--pipeline-script`). Files that appear by rename (rsync, `scp` to a temp
name, most sync tools) are detected as well as newly created ones. Runs never
overlap; files arriving during a run trigger one more run afterwards. Any
FASTQs already present at startup trigger one run immediately.

### Running in the background (nohup)

To keep the watcher running after you log out, start it with `nohup`. Run
the virtual environment's Python directly (not `uv run`), so the saved
process ID is the watcher itself and `kill` reaches it:

```bash
cd /mnt/volume/irma_pipe2
deactivate 2>/dev/null          # leave any other activated venv
uv sync                         # make sure .venv is up to date

nohup .venv/bin/python watch_and_run.py \
    -i /path/to/input \
    -o /path/to/results \
    --config config.json \
    --debounce 60 \
    > watch.log 2>&1 &

echo $! > watch.pid             # process ID, for stopping it later
```

Check on it:

```bash
tail -f watch.log                                  # live log (Ctrl+C stops viewing, not the watcher)
pgrep -af watch_and_run                            # is it running?
ls -t /mnt/volume/logs/pipeline_*.log | head -1    # newest pipeline run log (logs_dir)
```

Stop it:

```bash
kill $(cat watch.pid) && rm watch.pid
```

The watcher exits cleanly; a pipeline run that is already in progress
finishes on its own, so no sample is left half-done. Run only one watcher at
a time (check with `pgrep -af watch_and_run` before starting another).

A one-off run in the background, e.g. to reprocess a folder:

```bash
nohup .venv/bin/python run_pipeline.py -i /path/to/input -o /path/to/results \
    --config config.json > run_$(date +%Y%m%d_%H%M).log 2>&1 &
```

`nohup` does not restart the watcher after a server reboot or a crash; a
systemd service does both, if that is needed later.

## Output structure

```
<results>/<RUN_FOLDER>_RESULTS/
├── <sample>/
│   ├── fastqc_results/                   # Illumina
│   ├── <sample>_fastp_R1.fastq.gz, _R2.fastq.gz, _fastp_report.{html,json}   # Illumina
│   ├── <sample>.fastq.gz, nanoplot/      # MinION: merged reads, NanoPlot report
│   ├── kraken/<sample>_kraken_report.txt, _kraken_output.txt
│   ├── irma_results/<sample>/            # IRMA output: amended_consensus/, tables/, ...
│   └── snippy/                           # SARS-CoV-2, Illumina only
├── segments/<HA|NA|PB1|PB2|PA|MP|NS|NP>/<sample>_<SEGMENT>.fasta     # influenza
├── GISAIDfasta/<GISAID-name>/*.fasta, namelist.xlsx, submit_info.xlsx, fastas.xlsx   # influenza, Illumina
├── nextclade_results/<sample>/
├── nextclade_results.tsv                 # influenza: sample_name, Segmentti, Nextclade columns; one row per sample+segment
├── nextclade_results_COV.tsv             # SARS-CoV-2: sample_name + Nextclade columns
├── IRMA_analysis.xlsx                    # influenza N% per segment, read counts
├── sample_summary_report.xlsx            # Kraken, ambiguous sites, N%, coverage
└── multiqc_report.html, multiqc_data/    # refreshed whenever the run gets new results
```

### A note on IRMA's file names

IRMA names the files in `amended_consensus/` with its `SEG_NUMBERS`
setting (`IRMA_RES/scripts/amendConsensus.pl`): each `key:value` key is
matched against the gene name (`A_HA_H3`, `B_PB1`, …), the first match wins,
and the file becomes `<sample>_<value>.fa`. The lab uses

```
SEG_NUMBERS="B_PB1:B_PB1,B_PB2:B_PB2,A_PB2:A_PB2,A_PB1:A_PB1,PA:PA,HA:HA,NP:NP,NA:NA,MP:MP,NS:NS"
```

so influenza consensus files are `26-040001_HA.fa`, `26-040001_A_PB1.fa`,
`26-040001_B_PB2.fa`, … It is set in `IRMA_RES/modules/FLU/init.sh` **and**
in `custom-config-flu.sh`. The second copy matters: an IRMA upgrade replaces
`init.sh`, which would silently switch output back to IRMA's stock numbers.

`segments.py` turns these names into segments (HA, PB1, …) and every step
uses it. It also understands the stock numbering
(`B_PB1:1,B_PB2:2,A_PB2:1,A_PB1:2,PA:3,HA:4,NP:5,NA:6,M:7,NS:8`). There,
1 and 2 mean different genes for A and B, so the type is read from IRMA's
`A_*.fasta` / `B_*.fasta` files first.

Single-genome modules (CoV, RSV) have no `SEG_NUMBERS` and write
`<sample>.fa`. With `PADDED_CONSENSUS=1` (CoV and RSV defaults) IRMA also writes
an N-padded `<sample>.pad.fa`; the pipeline always ignores the padded copies.

### Nextclade output

The run-level files are tab-separated (`.tsv`): Nextclade's mutation columns
contain commas and semicolons, which broke the old CSV files. Read them with
a tab separator, e.g. `pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)`
(`keep_default_na=False` stops pandas turning the segment "NA" into a blank).
In Excel, open them via *Data → From Text/CSV* and choose Tab.

## Adding a species

Everything species-specific is declared in **`species.py`**. Reading the
input, QC, IRMA, archiving, done markers, the summary report and MultiQC are
shared. Example: RSV on MiSeq/NextSeq and MinION.

1. **Check IRMA has a module for it.** Look in `IRMA_RES/modules/` of your
   IRMA install ([CDCgov/irma](https://github.com/CDCgov/irma)). RSV has
   `RSV` and `config/minion.sh` → modules `RSV` and `RSV-minion`.
2. **Add an entry** to `SPECIES` in `species.py`:

   ```python
   from functools import partial
   from genome_pipeline import genome_post_irma

   SpeciesDef(
       code="RSV",                     # used in logs and config keys
       folder="rsv",                   # input/rsv/<sequencer>/<RUN_FOLDER>/
       irma_modules={"illumina": "RSV", "minion": "RSV-minion"},
       post_irma=partial(genome_post_irma,
                         nextclade_ref_key="nextclade_ref_rsv",
                         results_name="nextclade_results_RSV.tsv"),
       report_prefix="RSV",            # summary columns "RSV N%", "RSV average coverage", ...
   ),
   ```

   Leave a platform out of `irma_modules` if it isn't supported. Files for
   it are then left in the input folder, never archived unanalysed.
3. **Add its config keys** to `config.json`, here
   `"nextclade_ref_rsv": "/mnt/volume/db/nextclade/dataset/rsv_a"`. Add a
   `custom-config-rsv.sh` only if RSV needs IRMA settings of its own.
4. **Test it.** Copy an existing end-to-end test (or
   `tests/test_species.py::test_new_species_needs_only_a_registry_entry`),
   then run one real sample with `--dry-run` first.

`genome_post_irma` covers viruses with one genome and one Nextclade dataset.
A species that needs more writes its own post-IRMA function taking a
`workflow.SampleRun`, the way `covid_post_irma` adds Snippy and
`flu_post_irma` adds segments and GISAID naming. For RSV, the A/B dataset
choice would go there; `post_irma.get_dataset_prefix_for_sample` does the
same for influenza subtypes.

A new **sequencer** that behaves like an existing platform (e.g. iSeq →
Illumina) is one line in `SEQUENCER_PLATFORMS` in `species.py`.

## Development

```bash
uv run pytest          # unit tests + an end-to-end MiSeq smoke test
uv run ruff check .
```

`tests/test_miseq_end_to_end.py` runs the whole pipeline on a fake MiSeq
run folder. `tests/fake_tools/conda` stands in for `conda` and writes small
files shaped like real tool output, so the test checks how the steps fit
together (file names, folders, reports, archiving), not the biology. It
takes a few seconds and needs no conda environments.

If the end-to-end tests fail with "`conda_calls.log` not found", or pytest
stops with "is not executable", the fake `conda` has lost its execute
permission (common after uploading with WinSCP with "set permissions" off).
The tests would then call the real tools. Fix it with
`chmod +x tests/fake_tools/conda` (or `sudo chmod` if another user owns the
file); `tests/conftest.py` checks this before every test run.

Code layout:

| File | Role |
|---|---|
| `run_pipeline.py` | CLI, sample discovery, per-sample orchestration, archiving |
| `species.py` | **Species registry**: input folder, IRMA modules, post-IRMA step, report columns; sequencer → platform map |
| `workflow.py` | Shared per-sample workflow: reads → IRMA → species' post-IRMA step |
| `influenza_pipeline.py` | Influenza post-IRMA: segments, IRMA_analysis.xlsx, GISAID naming, Nextclade |
| `covid_pipeline.py` | SARS-CoV-2 post-IRMA: Snippy (Illumina), Nextclade |
| `genome_pipeline.py` | Generic post-IRMA for single-genome species (Nextclade) |
| `router.py` | Compatibility shim; routing now lives in `species.py` |
| `pre_irma.py` | FastQC / fastp / Kraken2 / Snippy / NanoPlot wrappers, MinION read merging, sample grouping, pairing checks |
| `irma.py` | IRMA wrapper |
| `post_irma.py` | Reports, Nextclade, MultiQC |
| `gisaid_naming.py`, `gisaid_naming_rerun.py` | GISAID naming (in-pipeline and standalone re-run) |
| `segments.py` | Segment names from IRMA file names |
| `watch_and_run.py` | Folder watcher |

## Known limitations

- **NextSeq runs with lane splitting** (`_L001`…`_L004` per sample): only the
  first R1/R2 pair is used (fastp gets two files); the other lanes are
  ignored. Run bcl2fastq/BCL Convert with `--no-lane-splitting`, or add a
  lane-merging step.
- **IRMA module names** in `species.py` (`FLU-minion`,
  `CoV-minion-long-reads`) are taken from IRMA's repository
  (`IRMA_RES/modules/*/config/`). If your install differs, override them with
  `irma_module_<code>_<platform>` in `config.json`.
- **MinION and GISAID**: MinION samples are never GISAID-named, by design;
  `gisaid_naming_rerun.py` also skips `barcode*` samples.
- **MinION reads are not filtered**: `fastq_fail` reads (lower quality) are
  used as well, by choice. Check IRMA's coverage and N% columns with that in mind.
- `nextclade_results.tsv` and `submit_info.xlsx` are appended to, so
  re-running a sample (`--force`) adds duplicate rows.

## Change log

### 2026-10-06 – MiSeq review

Every change is marked in the code with `CHANGE (2026-10-06)`.

Bugs that stopped influenza samples (NextSeq and MiSeq share the workflow):
- Segment copy crashed on every sample: `.pad.fa` twins (and IRMA's stock
  numbered names, `_4.fa`) were parsed as segments "HA.PAD" / "4".
  Because `.IRMA_DONE` had already been written, the sample was then
  skipped on all later runs, with no IRMA_analysis, GISAID name or Nextclade.
- `.IRMA_DONE` is now written only after the whole sample succeeded.
- N% columns, GISAID files and Nextclade now ignore `.pad.fa` copies and
  understand numbered names (N% was "No consensus" or inflated; GISAID
  submit columns were empty; Nextclade never found a dataset).
- Nextclade results are now TSV (`nextclade_results.tsv`, `nextclade_results_COV.tsv`) read from
  Nextclade's `nextclade.tsv`, with a `Segmentti` column; the COVID file also gets `sample_name`.
  The old CSV also lost the label of every NA-segment row (pandas reads "NA" as missing).

Raw-file handling:
- Raw files were archived even when a workflow stopped early (missing env,
  fastp or IRMA failure), for unsupported combinations (covid/miseq) and for
  the unimplemented MinION workflow. Now only completed samples are archived.
- `Undetermined` reads were analysed and GISAID-named as a sample; now skipped.
- `--dry-run` was ignored: it ran every tool and archived files. It now changes nothing.
- `archive_dir` from the config is used (was hard-coded); `--logs-dir` and
  `logs_dir` no longer disagree; IRMA's external config path is configurable.
- Running without `--config` exited silently; `--config` is now required.
- Results, per-sample logs, `.IRMA_DONE` markers and the archive are keyed
  by the sequencing-run folder instead of the date the pipeline started
  (`<RUN_FOLDER>_RESULTS`). A run processed on both sides of midnight (late
  R2, retry, watcher re-trigger) was split across two dated folders, with
  its reports, `segments/`, `GISAIDfasta/` and MultiQC split too.
- The same sample name in two different runs is no longer merged into one
  four-file "sample" (fastp used two of them and all four were archived).
- One log file per invocation; the day's `pipeline.log` used to be
  overwritten by each watcher-triggered run.
- MultiQC runs once per run folder that got new results.
- MiSeq COVID enabled (same workflow as NextSeq COVID), with an end-to-end test.
- MinION workflows for influenza and SARS-CoV-2: barcode reads (pass + fail)
  merged, NanoPlot, Kraken2 (single-end), IRMA `FLU-minion` / `CoV-minion-long-reads`,
  same post-processing as Illumina except GISAID naming (flu) and Snippy
  (CoV). `unclassified` reads skipped; pass/fail kept apart in the archive;
  R1/R2 pairing check applies to Illumina only; files are size-checked in
  one batch. End-to-end tests for both.
- `gisaid_naming_rerun.py` never names MinION `barcode*` samples (`--all`
  skips them, `--sample barcodeNN` refuses).
- Kraken2 now gets `--paired` for Illumina pairs (R1 and R2 were classified
  as unrelated single reads).
- **Species registry** (`species.py`) and one shared workflow (`workflow.py`)
  replace the four per-species/per-platform workflow functions, the
  routing table in `router.py`, the folder map in `pre_irma.py` and the
  `'FLU'`/`'CoV'` branches in the summary report. A new species is a
  registry entry; see "Adding a species". A test adds an RSV-like species
  this way and runs it end to end.
- Checked against IRMA's source (github.com/CDCgov/irma):
  - IRMA has no `CoV-minion` config, so the COVID MinION default is now
    `CoV-minion-long-reads`; the old default would have failed on the first run.
  - IRMA's stock segment numbers 1/2 differ for influenza A and B
    (`A_PB2:1`, `B_PB1:1`); `segments.py` now decodes them by type instead of
    assuming A.
  - The lab's `SEG_NUMBERS` is now also in `custom-config-flu.sh`, an
    influenza-only external config, so IRMA upgrades can't reset it and it
    doesn't apply to other modules.

Other:
- `openpyxl` added to dependencies (needed by `pandas.read_excel`); `ruff` moved to dev.
- R1/R2 detection matches only the end of the file name.
- Watcher reacts to renamed-in files (rsync etc.), not only newly created ones.
- A missing Kraken2 env/db no longer skips IRMA; a GISAID-naming error no
  longer skips Nextclade; a missing Snippy env no longer skips Nextclade.
- Snippy's arguments were passed in the wrong order (worked by accident).
- MultiQC runs with `--force` (no more `multiqc_report_1.html`, `_2`, …).
- Simplified: shared `run_pre_irma_steps()` (was duplicated in both
  workflows), `segments.py` (replaces three different parsers),
  `gisaid_naming.Sample.get_segments` (dead fallback branch removed),
  `segments_to_submit` (if/elif chain → lookup), `run_pipeline.main()` split
  into `parse_args`, `setup_logging`, `process_sample`, `archive_sample_files`.
- New tests in `tests/`, including the end-to-end MiSeq test.
- `tests/conftest.py` (2026-10-15) stops the test run if the fake `conda`
  is not executable, instead of silently running the real tools.
- README: running the watcher in the background with `nohup` (2026-10-15).

### 2026-09-23
- Samples processed in parallel (`--max-workers`).
- Incomplete or still-growing R1/R2 pairs are left in place instead of archived.
- One sample's crash no longer stops the run; empty segment consensus no longer crashes.

## Suggested next refactors

1. `gisaid_naming._process_sample_locked` (~300 lines) still mixes type and
   subtype detection, three report files and the central list update, with
   several broad `except Exception` blocks. Split into
   `detect_type_and_subtype()`, `write_namelist()`, `write_submit_info()`,
   `write_fastas_list()`, `update_central_list()`, and let real errors show.
2. `run_nexclade_analysis` and `run_nextclade_covid_analysis` share most of
   their code; a common `_run_nextclade_and_merge()` would remove the duplication.
3. `post_irma.Sample` stores results in a 17-item list addressed by index
   (`result[14]`); a dict keyed by `HEADER` names would be safer.
4. Validate the config once at startup (required envs exist, paths exist)
   instead of discovering missing entries sample by sample.