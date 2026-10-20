# Outputs and logs

[← README](../README.md)

## Results folder

Results go to `<results>/RESULTS_<YYYYMMDD>/`, named after the date the
pipeline invocation started. Every invocation on the same day (watcher
re-runs, retries) adds to the same folder and its run-level reports. A sample
completed by an invocation after midnight, e.g. because its R2 arrived late,
goes to the next day's folder.

```
<results>/RESULTS_<YYYYMMDD>/
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

`nextclade_results.tsv` and `submit_info.xlsx` are appended to, so re-running
a sample with `--force` adds duplicate rows.

## Logs, done markers and archive

| What | Where |
|---|---|
| One log per invocation | `<logs_dir>/pipeline_<YYYYMMDD_HHMMSS>.log` |
| Per-sample log with full tool output | `<logs_dir>/RESULTS_<YYYYMMDD>/<sample>/pipeline.log` |
| Done marker | `<logs_dir>/RESULTS_<YYYYMMDD>/<sample>/.IRMA_DONE` |
| Archived raw FASTQs | `<archive_dir>/RESULTS_<YYYYMMDD>/<sample>/` (MinION: `…/<sample>/fastq_pass/` and `fastq_fail/`) |

A sample with a `.IRMA_DONE` marker is skipped on later runs **the same day**
unless `--force` is given. The marker lives in that day's folder, so if the
same sample's FASTQs are copied in again on a later day, it is processed
again into that day's results.

## Reading the Nextclade TSVs

The run-level files are tab-separated (`.tsv`): Nextclade's mutation columns
contain commas and semicolons, which broke the old CSV files. Read them with
a tab separator:

```python
pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
```

`keep_default_na=False` stops pandas turning the segment "NA" into a blank.
In Excel, open them via *Data → From Text/CSV* and choose Tab.