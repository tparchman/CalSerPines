### Notes on cleaning NovaSeq flow cell from 11/20, barcode parsing for multiple projects, and analyses of Attenuata data 


`NOTE`: Here we need to clean two lanes of NovaSeq data. The first (sample-1_S1_L001_R1_001.fastq.gz) has PHHA, ACTH, CHDO, and PIAT. The second (sample-2_S2_L002_R1_001.fastq.gz) has POSE, PIAT, and some quadrus. We will need to clean and parse both to get all of the attenuata samples into working shape. 

`NOTE`: Most analyses will be completed at 

### This file contains code and notes for
1) cleaning contaminants using tapioca
2) parsing barcodes
3) splitting fastqs 
4) 
5) calling variants
6) filtering
7) entropy for genotype probabilities.

## Cleaning contaminants

Being executed on ponderosa using tapioca pipeline. Commands in bash script, executed as below (11/17/20).

    $ module load fqutils/0.4.1
    $ module load bowtie2/2.2.5
    
    $ bash cleaning_bash.sh &


After .clean.fastq has been produced, recompress raw data:

    $ gzip sample-1_S1_L001_R1_001.fastq &
    $ gzip sample-2_S2_L002_R1_001.fastq &


Number of reads **before** cleaning:

    $ grep -c "^@" sample-1_S1_L001_R1_001.fastq > S1_number_of_rawreads.txt
    $ grep -c "^@" sample-2_S2_L002_R1_001.fastq > S2_number_of_rawreads.txt
    
Number of reads **after** cleaning:

    $ grep "^@" S1_11_20.clean.fastq -c > S1_No_ofcleanreads.txt &
    $ less S1_No_ofcleanreads.txt
    # lane 1: 1,520,946,461

    $ grep "^@" S2_11_20.clean.fastq -c > S2_No_ofcleanreads.txt &
    $ less S2_No_ofcleanreads.txt
    # lane 2: 

# DONE TO HERE, BELOW IS PLACEHOLDER FROM TPOD

## Barcode parsing:

Barcode keyfile are `/mnt/UTGSAF_11_20/11_20_GSAF_lane1BCODEKEY.csv` and `/mnt/UTGSAF_11_20/11_20_GSAF_lane1BCODEKEY.csv`

Parsing library 1:

    $ perl parse_barcodes768.pl 11_20_GSAF_lane1BCODEKEY.csv S1_11_20.clean.fastq A00 &

# DONE TO HERE, BELOW IS PLACEHOLDER FROM TPOD

Parsing library 2:

    $ perl parse_barcodes768.pl 11_20_GSAF_lane2BCODEKEY.csv S2_11_20.clean.fastq A00 &

`NOTE`: the A00 object is the code that identifies the sequencer (first three characters after the @ in the fastq identifier).

    $ less parsereport_tpodura.clean.fastq
    #Good mids count: 1617562664
    #Bad mids count: 58132389
    #Number of seqs with potential MSE adapter in seq: 305112
    #Seqs that were too short after removing MSE and beyond: 193
          
Cleaning up the directory:

    $ rm tpodura.clean.fastq
    $ rm miderrors_tpodura.clean.fastq
    $ rm parsereport_tpodura.clean.fastq

### Alignment to *T. cristinae* genome and variant calling.
New software needs to be installed on ponderosa, including:
- bwa 0.7.17-r1188 (https://github.com/lh3/bwa/releases)
- bcftools 1.9 (under https://sourceforge.net/projects/samtools/files/samtools/1.9/)
- samtools 1.10 (under https://sourceforge.net/projects/samtools/files/samtools/1.10/)
## mapping with `bwa`


## alignment with bwa 0.7.17-r1188
## assumes genome = genome.fasta and the sample is tpod1

bwa mem -k 20 -w 100 -r 1.3 -T 30 -R @RG        ID:tpod1        PL:ILLUMINA     LB:tpod1        SM:tpod1 genome.fasta tpod1.fastq > aln_tpod1.sam 2> tpod1.log

## standard sam to bam, sorting and indexing; code not included, same as it has always been
## at present using samtools 1.10 and bcftools 1.9

## variant calling with bcftools 1.9
## bams is a text file with all of the sorted bam files listed, one per line
bcftools mpileup -C 50 -d 250 -f genome.fasta -q 30 -Q 20 -I -b bams -O b -o tpod.bcf

## sometimes I use the -c option, sometimes not, I have mixed feelings, probably would use it
## for podura
bcftools call -v -c -f GQ -p 0.01 -P 0.001 -O v -o tpod.vcf tpod.bcf
