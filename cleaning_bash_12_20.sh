#!/bin/bash

##########################################################################################
## UTGSAF 12_20 cleaning 1 (bighorn sheep,PIMU)
##########################################################################################

/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/illumina_oligos --pct 20 /archive/parchman_lab/rawdata_to_backup/GSAF_12_20/BH-sheep-ddRADseq_S1_L001_R1_001.fastq > BHS_PIMU.readstofilter.ill.txt 

echo "Illumina filtering done for lane 1"

/working/jahner/tapioca/src/tap_contam_analysis --db /archive/parchman_lab/rawdata_to_backup/contaminants/phix174 --pct 80 /archive/parchman_lab/rawdata_to_backup/GSAF_12_20/BH-sheep-ddRADseq_S1_L001_R1_001.fastq  > BHS_PIMU.readstofilter.phix.txt 

echo "PhiX filtering done for lane 1"


/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/ecoli-k-12 --pct 80 /archive/parchman_lab/rawdata_to_backup/GSAF_12_20/BH-sheep-ddRADseq_S1_L001_R1_001.fastq  > BHS_PIMU.readstofilter.ecoli.txt

echo "ecoli filtering done for lane 1"


cat /archive/parchman_lab/rawdata_to_backup/GSAF_12_20/BH-sheep-ddRADseq_S1_L001_R1_001.fastq | fqu_cull -r BHS_PIMU.readstofilter.ill.txt BHS_PIMU.readstofilter.phix.txt BHS_PIMU.readstofilter.ecoli.txt > BHS_PIMU.clean.fastq

echo "Clean copy of lane 1 done"

##########################################################################################
## UTGSAF 12_20 cleaning 2 (TICR and PIMU)
##########################################################################################

/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/illumina_oligos --pct 20 /archive/parchman_lab/rawdata_to_backup/GSAF_12_20/Tcristinae-ddRADseq_S2_L002_R1_001.fastq > TICR_PIMU.readstofilter.ill.txt 

echo "Illumina filtering done for lane 2"

/working/jahner/tapioca/src/tap_contam_analysis --db /archive/parchman_lab/rawdata_to_backup/contaminants/phix174 --pct 80 /archive/parchman_lab/rawdata_to_backup/GSAF_12_20/Tcristinae-ddRADseq_S2_L002_R1_001.fastq  > TICR_PIMU.readstofilter.phix.txt 

echo "PhiX filtering done for lane 2"


/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/ecoli-k-12 --pct 80 /archive/parchman_lab/rawdata_to_backup/GSAF_12_20/Tcristinae-ddRADseq_S2_L002_R1_001.fastq  > TICR_PIMU.readstofilter.ecoli.txt

echo "ecoli filtering done for lane 2"


cat /archive/parchman_lab/rawdata_to_backup/GSAF_12_20/Tcristinae-ddRADseq_S2_L002_R1_001.fastq | fqu_cull -r TICR_PIMU.readstofilter.ill.txt TICR_PIMU.readstofilter.phix.txt TICR_PIMU.readstofilter.ecoli.txt > TICR_PIMU.clean.fastq

echo "Clean copy of lane 2 done"


