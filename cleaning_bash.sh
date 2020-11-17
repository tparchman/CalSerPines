#!/bin/bash

##########################################################################################
## T.podura cleaning 1
##########################################################################################

/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/illumina_oligos --pct 20 /working/parchman/tpodura/L2-timema_S2_L002_R1_001.fastq > tpodura.readstofilter.ill.txt 

echo "Illumina filtering done for lane 2"

/working/jahner/tapioca/src/tap_contam_analysis --db /archive/parchman_lab/rawdata_to_backup/contaminants/phix174 --pct 80 /working/parchman/tpodura/L2-timema_S2_L002_R1_001.fastq > tpodura.readstofilter.phix.txt 

echo "PhiX filtering done for lane 2"


/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/ecoli-k-12 --pct 80 /working/parchman/tpodura/L2-timema_S2_L002_R1_001.fastq > tpodura.readstofilter.ecoli.txt

echo "ecoli filtering done for lane 2"


cat /working/parchman/tpodura/L2-timema_S2_L002_R1_001.fastq | fqu_cull -r tpodura.readstofilter.ill.txt tpodura.readstofilter.phix.txt tpodura.readstofilter.ecoli.txt > tpodura.clean.fastq

echo "Clean copy of lane 2 done"


