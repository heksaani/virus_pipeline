# IRMA consensus file names

[← README](../README.md)

IRMA names the files in `amended_consensus/` with its `SEG_NUMBERS`
setting (`IRMA_RES/scripts/amendConsensus.pl`): each `key:value` key is
matched against the gene name (`A_HA_H3`, `B_PB1`, …), the first match wins,
and the file becomes `<sample>_<value>.fa`.

## The lab's setting

```
SEG_NUMBERS="B_PB1:B_PB1,B_PB2:B_PB2,A_PB2:A_PB2,A_PB1:A_PB1,PA:PA,HA:HA,NP:NP,NA:NA,MP:MP,NS:NS"
```

so influenza consensus files are `26-040001_HA.fa`, `26-040001_A_PB1.fa`,
`26-040001_B_PB2.fa`, …

It is set in `IRMA_RES/modules/FLU/init.sh` **and** in
`custom-config-flu.sh`. The second copy matters: an IRMA upgrade replaces
`init.sh`, which would silently switch output back to IRMA's stock numbers.
`tests/test_species.py::test_flu_external_config_keeps_lab_seg_numbers`
fails if the line is removed from `custom-config-flu.sh`.

## IRMA's stock numbering

`segments.py` turns these names into segments (HA, PB1, …) and every step
uses it. It also understands the stock numbering:

```
B_PB1:1,B_PB2:2,A_PB2:1,A_PB1:2,PA:3,HA:4,NP:5,NA:6,M:7,NS:8
```

There, 1 and 2 mean different genes for influenza A and B, so the type is
read from IRMA's `A_*.fasta` / `B_*.fasta` files first.

## Single-genome modules

Single-genome modules (CoV, RSV) have no `SEG_NUMBERS` and write
`<sample>.fa`. With `PADDED_CONSENSUS=1` (CoV and RSV defaults) IRMA also
writes an N-padded `<sample>.pad.fa`; the pipeline always ignores the padded
copies.