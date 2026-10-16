# IRMA external config for influenza (FLU and FLU-minion).
# Used instead of custom-config.sh for influenza samples (see irma.irma_external_config).
# Keep the shared settings below in step with custom-config.sh.
MIN_CONS_SUPPORT=20
TMP=/mnt/volume/irma_pipe/temp_irma/

# Name amended_consensus files by segment (<sample>_HA.fa, <sample>_A_PB1.fa, ...)
# instead of IRMA's stock numbers (<sample>_4.fa). This is the lab's setting from
# flu-amd/IRMA_RES/modules/FLU/init.sh; set here too so an IRMA upgrade, which
# replaces init.sh, cannot silently switch the output back to numbers.
SEG_NUMBERS="B_PB1:B_PB1,B_PB2:B_PB2,A_PB2:A_PB2,A_PB1:A_PB1,PA:PA,HA:HA,NP:NP,NA:NA,MP:MP,NS:NS"
