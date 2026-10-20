# Changelog

## 2026-10-20 – Documentation split

- README cut down to setup and everyday use. Reference material moved to
  `docs/` (configuration, input layout, outputs, IRMA file naming, adding a
  species), test and code-layout notes to `CONTRIBUTING.md`, this change log
  to `CHANGELOG.md`. Running the watcher on the server and design rationale
  moved to the Gitea wiki; the "Suggested next refactors" list became Gitea
  issues.

## 2026-10-15

- `tests/conftest.py` stops the test run if the fake `conda` is not
  executable, instead of silently running the real tools.
- README: running the watcher in the background with `nohup`.

## 2026-10-06 – MiSeq review

Every change is marked in the code with `CHANGE (2026-10-06)`.

After pulling this change, run `uv lock` and `uv sync` (`openpyxl` was added
to `pyproject.toml`).

### Bugs that stopped influenza samples (NextSeq and MiSeq share the workflow)

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

### Raw-file handling

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

### New platforms and species structure

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
  registry entry; see `docs/adding-a-species.md`. A test adds an RSV-like
  species this way and runs it end to end.
- Checked against IRMA's source (github.com/CDCgov/irma):
  - IRMA has no `CoV-minion` config, so the COVID MinION default is now
    `CoV-minion-long-reads`; the old default would have failed on the first run.
  - IRMA's stock segment numbers 1/2 differ for influenza A and B
    (`A_PB2:1`, `B_PB1:1`); `segments.py` now decodes them by type instead of
    assuming A.
  - The lab's `SEG_NUMBERS` is now also in `custom-config-flu.sh`, an
    influenza-only external config, so IRMA upgrades can't reset it and it
    doesn't apply to other modules.

### Other

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

## 2026-09-23

- Samples processed in parallel (`--max-workers`).
- Incomplete or still-growing R1/R2 pairs are left in place instead of archived.
- One sample's crash no longer stops the run; empty segment consensus no longer crashes.