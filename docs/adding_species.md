# Adding a species or sequencer

[← README](../README.md)

Everything species-specific is declared in **`species.py`**. Reading the
input, QC, IRMA, archiving, done markers, the summary report and MultiQC are
shared. Example: RSV on MiSeq/NextSeq and MinION.

## Steps

1. **Check IRMA has a module for it.** Look in `IRMA_RES/modules/` of your
   IRMA install ([CDCgov/irma](https://github.com/CDCgov/irma)). RSV has
   `RSV` and `config/minion.sh` → modules `RSV` and `RSV-minion`.
2. **Add an entry** to `SPECIES` in `species.py`:

   ```python
   from functools import partial
   from genome_pipeline import genome_post_irma

   SpeciesDef(
       code="RSV",                     # used in logs and config keys
       folder="rsv",                   # input/rsv/<sequencer>/...
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

## Species that need more than Nextclade

`genome_post_irma` covers viruses with one genome and one Nextclade dataset.
A species that needs more writes its own post-IRMA function taking a
`workflow.SampleRun`, the way `covid_post_irma` adds Snippy and
`flu_post_irma` adds segments and GISAID naming. For RSV, the A/B dataset
choice would go there; `post_irma.get_dataset_prefix_for_sample` does the
same for influenza subtypes.

## Adding a sequencer

A new sequencer that behaves like an existing platform (e.g. iSeq →
Illumina) is one line in `SEQUENCER_PLATFORMS` in `species.py`.