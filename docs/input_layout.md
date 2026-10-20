# Input layout

[← README](../README.md)

```
<input_dir>/<influenza|covid>/<nextseq|miseq|minion>/.../<sample>_R1.fastq.gz
```

The first folder under the input folder is the species and the second is the
sequencer; both are case-insensitive. Below them the pipeline searches every
subfolder, so any layout works.

## Illumina (NextSeq, MiSeq)

Copy one folder per sample, or just the FASTQs:

```
input/influenza/miseq/
    26-040001/
        26-040001_R1.fastq.gz
        26-040001_R2.fastq.gz
    26-040002_S2_L001_R1_001.fastq.gz         <- loose FASTQs work too
    26-040002_S2_L001_R2_001.fastq.gz
    Undetermined_S0_L001_R1_001.fastq.gz      <- skipped, left in place
```

A whole MiSeq run folder (`.../Data/Intensities/BaseCalls/*.fastq.gz`) also
works, but the run folder name is not used for anything.

When a sample completes, its FASTQs are archived and a sample folder left
empty is removed.

### File names

- **R1/R2** are read from the end of the file name: `_R1_001.fastq.gz`,
  `_R1.fastq.gz`, `_1.fq.gz` and the same with 2.
- **Extensions:** `.fastq.gz`, `.fq.gz`, `.fastq` or `.fq`. A file named
  `26-040001_R1.gz` is not found at all.

### Sample names

The sample name comes from the **file name**, not the folder: the part
before `_S<number>` (`26-040002_S2_L001_R1_001.fastq.gz` → `26-040002`), or
if there is none, the part before `_R1`/`_R2` (`26-040001_R1.fastq.gz` →
`26-040001`). It is also the THL sample name used for GISAID naming, so it
must be in `YY-NNNNNN` or `YYYY-NNNN` form. `Undetermined` reads are never
analysed.

> ⚠️ **Each sample name only once.** All files with the same sample name are
> treated as one sample, wherever they are under the input folder, even in a
> different species or sequencer folder. The species and sequencer are then
> taken from one of the files, and fastp uses only the first two. Don't copy
> a re-sequenced sample in while the earlier copy is still in the input
> folder.

### Lane splitting

If a NextSeq run is demultiplexed with lane splitting (`_L001`…`_L004` per
sample), only the first R1/R2 pair is used and the other lanes are ignored.
Run bcl2fastq / BCL Convert with `--no-lane-splitting`, or add a
lane-merging step.

## MinION

Copy the finished run folder (the one MinKNOW writes, with `fastq_pass/` and
`fastq_fail/`) under `<species>/minion/`:

```
input/influenza/minion/20261006_1203_MN41234_FAX12345_abcd1234/
    fastq_pass/barcode01/*.fastq.gz      <- sample "barcode01"
    fastq_fail/barcode01/*.fastq.gz      <- same sample, also used
    fastq_pass/unclassified/             <- skipped, left in place
```

Each `barcodeNN` folder is one sample, named `barcodeNN` throughout the
results. Gzipped and plain FASTQs are both accepted. Because these are not
THL sample names, MinION samples get no GISAID name.

Barcode names repeat between runs, so process one MinION run at a time: two
runs' `barcode01` folders in the input at once are merged into one sample.

The pipeline does not wait for MinKNOW to finish, so **only copy runs that
have completed**.