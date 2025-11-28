FROM cdcgov/irma

# Replace the numbered segments with named segments
RUN sed -i 's/SEG_NUMBERS="B_PB1:1,B_PB2:2,A_PB2:1,A_PB1:2,PA:3,HA:4,NP:5,NA:6,M:7,NS:8"/SEG_NUMBERS="B_PB1:B_PB1,B_PB2:B_PB2,A_PB2:A_PB2,A_PB1:A_PB1,PA:PA,HA:HA,NP:NP,NA:NA,M:MP,NS:NS"/' \
    /app/IRMA_RES/modules/FLU/init.sh