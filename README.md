# Virus pipeline

Pipeline for viral sequencing data


## IRMA

Use the dockerfile to create container images

``` bash
docker build -t irma_thl .
```

Use the image:
``` bash
docker run --rm -it \
    -v $(pwd):/data \
    ghcr.io/cdcgov/irma:latest \
    IRMA FLU sample_R1_001.fastq.gz sample_R2_001.fastq.gz ./sample1
```