# Virus pipeline

Pipeline for viral sequencing data

```mermaid
  graph TD;
      minION--> barcode.fastq;
      nextSeq-->R1.fastq1,R2.fastq;
      Miseq-->R1.fastq1,R2.fastq;
      barcode.fastq-->combine;
      combine-->IRMA;
      R1.fastq1,R2.fastq-->IRMA;
      IRMA-->seqments_to_folders.py;
      seqments_to_folders.py-->IRMA_analysis_with_reads.py;
      IRMA_analysis_with_reads.py-->Infl_GISAID_naming.py
```

## IRMA

Use the dockerfile to create container image

``` bash
docker build -t irma_thl .
```

Use the image to analyze the sample 
``` bash
docker run --rm -it \
    -v $(pwd):/data \
    irma_thl \
    IRMA FLU sample_R1_001.fastq.gz sample_R2_001.fastq.gz ./sample1
```
