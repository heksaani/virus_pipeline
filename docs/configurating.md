# Configuration

[← README](../README.md)

The pipeline reads one JSON file, given with `--config`.

## `config.json` keys

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

## IRMA external configs

IRMA reference: <https://github.com/CDCgov/irma> (modules and their configs
are in `IRMA_RES/modules/<MODULE>/`).

`custom-config.sh` sets IRMA options (e.g. `MIN_CONS_SUPPORT=20`) and IRMA's
`TMP` folder. Check that `TMP` points at a folder that exists on your server.
Influenza uses `custom-config-flu.sh` instead: the same settings plus the
lab's `SEG_NUMBERS` (see [irma-naming.md](irma-naming.md)). Keep the shared
lines of the two files in step. IRMA applies the external config after the
module's own `init.sh` and config, so these values win.

## Nextclade datasets

`nextclade_ref` must contain one folder per subtype and segment, named
`flu_<prefix>_<segment>`, e.g. `flu_h3n2_ha`, `flu_h1n1pdm_na`, `flu_b_ha`.
Segments without a matching folder are skipped with a warning.

`nextclade_ref_covid` points directly at the SARS-CoV-2 dataset folder.