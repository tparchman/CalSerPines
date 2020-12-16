#!/bin/bash

##########################################################################################
## UTGSAF 11_20 cleaning 1 (PHHA, ACTH, CHDO, PIAT)
##########################################################################################

/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/illumina_oligos --pct 20 /mnt/UTGSAF_11_20/sample-1_S1_L001_R1_001.fastq > S1_11_20.readstofilter.ill.txt 

echo "Illumina filtering done for lane 1"

/working/jahner/tapioca/src/tap_contam_analysis --db /archive/parchman_lab/rawdata_to_backup/contaminants/phix174 --pct 80 /mnt/UTGSAF_11_20/sample-1_S1_L001_R1_001.fastq  > S1_11_20.readstofilter.phix.txt 

echo "PhiX filtering done for lane 1"


/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/ecoli-k-12 --pct 80 /mnt/UTGSAF_11_20/sample-1_S1_L001_R1_001.fastq  > S1_11_20.readstofilter.ecoli.txt

echo "ecoli filtering done for lane 1"


cat /mnt/UTGSAF_11_20/sample-1_S1_L001_R1_001.fastq | fqu_cull -r S1_11_20.readstofilter.ill.txt S1_11_20.readstofilter.phix.txt S1_11_20.readstofilter.ecoli.txt > S1_11_20.clean.fastq

echo "Clean copy of lane 1 done"

##########################################################################################
## UTGSAF 11_20 cleaning 2 (POSE, Quadrus, PIAT)
##########################################################################################

/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/illumina_oligos --pct 20 /mnt/UTGSAF_11_20/sample-2_S2_L002_R1_001.fastq > S2_11_20.readstofilter.ill.txt 

echo "Illumina filtering done for lane 2"

/working/jahner/tapioca/src/tap_contam_analysis --db /archive/parchman_lab/rawdata_to_backup/contaminants/phix174 --pct 80 /mnt/UTGSAF_11_20/sample-2_S2_L002_R1_001.fastq  > S2_11_20.readstofilter.phix.txt 

echo "PhiX filtering done for lane 2"


/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/ecoli-k-12 --pct 80 /mnt/UTGSAF_11_20/sample-2_S2_L002_R1_001.fastq  > S2_11_20.readstofilter.ecoli.txt

echo "ecoli filtering done for lane 2"


cat /mnt/UTGSAF_11_20/sample-2_S2_L002_R1_001.fastq | fqu_cull -r S2_11_20.readstofilter.ill.txt S2_11_20.readstofilter.phix.txt S2_11_20.readstofilter.ecoli.txt > S2_11_20.clean.fastq

echo "Clean copy of lane 2 done"


