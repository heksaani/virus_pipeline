# Contributing

## Setup

```bash
uv sync          # installs pytest and ruff as dev dependencies
```

## Tests and linting

```bash
uv run pytest          # unit tests + end-to-end smoke tests
uv run ruff check .
```

The end-to-end tests (`tests/test_miseq_end_to_end.py`,
`tests/test_covid_miseq_end_to_end.py`, `tests/test_minion_end_to_end.py`) run
the whole pipeline on fake run folders. `tests/fake_tools/conda` stands in
for `conda` and writes small files shaped like real tool output, so the tests
check how the steps fit together (file names, folders, reports, archiving),
not the biology. They take a few seconds and need no conda environments.

### "conda_calls.log not found" or "is not executable"

The fake `conda` has lost its execute permission (common after uploading with
WinSCP with "set permissions" off). The tests would then call the real tools.
Fix it with:

```bash
chmod +x tests/fake_tools/conda      # or sudo chmod if another user owns the file
```

`tests/conftest.py` checks this before every test run.

## Code layout

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

To add a virus or sequencer, see [docs/adding-a-species.md](docs/adding-a-species.md).

## Conventions

- Mark changes in the code with `CHANGE (YYYY-MM-DD)` and add an entry to
  [CHANGELOG.md](CHANGELOG.md).
- Keep documentation that describes code behaviour in `docs/`, in the same
  commit as the code change. Operational notes and design rationale go in
  the Gitea wiki.
- Planned refactors are Gitea issues, not README sections.