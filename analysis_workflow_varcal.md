### Notes on cleaning NovaSeq flow cell from 11/20, barcode parsing for multiple projects, and analyses of pine data 


`NOTE`: Here we need to clean two lanes of NovaSeq data. The first (sample-1_S1_L001_R1_001.fastq.gz) has PHHA, ACTH, CHDO, and PIAT. The second (sample-2_S2_L002_R1_001.fastq.gz) has POSE, PIAT, and some quadrus. We will need to clean and parse both to get all of the attenuata samples into working shape. 

`NOTE`: Most analyses will be completed at 

### This file contains code and notes for
1) cleaning contaminants using tapioca
2) parsing barcodes
3) splitting fastqs 
4) de novo assembly
5) reference based assembly
6) calling variants
7) filtering
8) entropy for genotype probabilities.

####################################################################################
## 1. Cleaning contaminants
####################################################################################

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
    # lane 2: 1,521,453,849

####################################################################################
## 2. Barcode parsing:
####################################################################################

Barcode keyfile are `/mnt/UTGSAF_11_20/11_20_GSAF_lane1BCODEKEY.csv` and `/mnt/UTGSAF_11_20/11_20_GSAF_lane1BCODEKEY.csv`

Parsing library 1:

    $ perl parse_barcodes768.pl 11_20_GSAF_lane1BCODEKEY.csv S1_11_20.clean.fastq A00 &

Parsing library 2:

    $ perl parse_barcodes768.pl 11_20_GSAF_lane2BCODEKEY.csv S2_11_20.clean.fastq A00 &

`NOTE`: the A00 object is the code that identifies the sequencer (first three characters after the @ in the fastq identifier).

    $ less parsereport_S1_11_20.clean.fastq
    Good mids count: 1484409965
    Bad mids count: 36536312
    Number of seqs with potential MSE adapter in seq: 337486
    Seqs that were too short after removing MSE and beyond: 184

    $ less parsereport_S1_11_20.clean.fastq
    Good mids count: 1475831187
    Bad mids count: 45622449
    Number of seqs with potential MSE adapter in seq: 332950
    Seqs that were too short after removing MSE and beyond: 213


Cleaning up the LIB1 directory:

    $ rm S1_11_20.clean.fastq
    $ rm miderrors_S1_11_20.clean.fastq
    $ rm parsereport_S1_11_20.clean.fastq

Cleaning up the LIB2 directory:

    $ rm S2_11_20.clean.fastq
    $ rm miderrors_S2_11_20.clean.fastq
    $ rm parsereport_S2_11_20.clean.fastq
    
####################################################################################
## 3. splitting fastqs
####################################################################################

Make ids file

    $ cut -f 3 -d "," 11_20_GSAF_lane1BCODEKEY.csv | grep "_" > L1_ids_noheader.txt

    $ cut -f 3 -d "," 11_20_GSAF_lane2BCODEKEY.csv | grep "_" > L2_ids_noheader.txt

Split fastqs by individual

    $ perl splitFastq_universal_regex.pl L1_ids_noheader.txt parsed_S1_11_20.clean.fastq &

    $ perl splitFastq_universal_regex.pl L2_ids_noheader.txt parsed_S2_11_20.clean.fastq &


Zip the parsed*fastq files for now, but delete once patterns and qc are verified:

    $ gzip parsed_S1_11_20.clean.fastq
    $ gzip parsed_S2_11_20.clean.fastq

Total reads for muricata, radiata, and attenuata (473 individuals)

    $ grep -c "^@" *fastq > seqs_per_ind.txt

### Moving fastqs to project specific directories

For LIB1:

    $mv AT_*fastq ACTH
    $mv PH_*fastq PHHA
    $mv PA_*fastq CalSer
    $mv CD_*fastq CHDO

For LIB2:

    $mv PS_*fastq poa
    $mv PX_*fastq CalSer
    $mv PA_*fastq CalSer
    $mv PM_*fastq CalSer
    $mv PR_*fastq CalSer
    $mv ONLY QUADRUS SHOULD BE LEFT.

### Moved the parsed files for all of the above project to /archive/parchman_lab/rawdata_to_backup

####################################################################################
### Lanie started here: begin with attenuata
####################################################################################

    /working/lgalland/attenuata/

####################################################################################
## 4. de novo assembly
####################################################################################


### removed individuals with <35M of sequence data
kept data threshold low to retain Mexican island populations

    rm -rf PA_LA_0074.fastq.gz
    rm -rf PA_OC_0093.fastq.gz
    rm -rf PA_OC_0096.fastq.gz
    rm -rf PA_OC_0100.fastq.gz
    rm -rf PA_SB_0082.fastq.gz
    rm -rf PA_SB_0083.fastq.gz
    rm -rf PA_SB_0085.fastq.gz
    rm -rf PM_DC_0236.fastq.gz
    rm -rf PM_FR_0019.fastq.gz
    rm -rf PM_NR_0064.fastq.gz
    rm -rf PM_SP_0210.fastq.gz
    rm -rf PM_SR_0037.fastq.gz
    rm -rf PR_CN_0819.fastq.gz
    rm -rf PR_CS_0822.fastq.gz
    rm -rf PR_GU_0818.fastq.gz
    rm -rf PR_GU_0820.fastq.gz
    rm -rf PR_GU_0821.fastq.gz
    rm -rf PR_MR_0780.fastq.gz
    rm -rf PX_SC_0810.fastq.gz

removed 19 individuals. new total = 454 individuals

### Combined pines, de novo assembly

#### 1. gzip files (already done, takes time)
    nohup gzip *fastq &>/dev/null &

#### 2. make list of fastqs
    ls *.fastq.gz > namelist

#### 3. one-liner to remove '.fastq.gz' from namelist (instant)
    sed -i'' -e 's/.fastq.gz//g' namelist
    
#### 4. set variables (run as one chunk command, instant)
    AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
    AWK2='!/>/'
    AWK3='!/NNN/'
    PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
    
#### 5. create list of sequences and counts for each individual (few minutes)
    cat namelist | parallel --no-notice -j 8 "zcat {}.fastq | mawk '$AWK1' | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs" &
    
#### 6. combine files for all individuals (about 15 seconds)
    cat *.uniq.seqs > uniq.seqs

#### 7. select those sequences that have 4 reads (or however you want to choose - pretty instant, about 1 minute)
    parallel --no-notice -j 8 mawk -v x=4 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp;    $z{$_}++;} while(($k,$v) = each(%z))     {print "$v\t$k\n";}' > uniqCperindv
    
    wc -l uniqCperindv
12696842 uniq sequences per individual

#### 8. restrict the data by how many individuals have that read
    for ((i = 2; i <= 10; i++));
    do
    echo $i >> ufile
    done

NOTE: the above bash loop just makes a list from 2 to 10
    
    cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data
    
    more uniqseq.peri.data

2	3719956
3	2020056
4	1347546
5	996134
6	780431
7	635769
8	532158
9	453801
10	392907

    rm -rf ufile
    
#### 9. restrict sequences to those that are found in at least 4 inds (few seconds per ind)
    mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.4.seqs
    wc -l uniq.k.4.c.4.seqs
1347546 uniq.k.4.c.4.seqs

#### 10. Convert these sequences to fasta format (both pretty instant)
    cut -f2 uniq.k.4.c.4.seqs > totaluniqseq
    mawk '{c= c + 1; print ">Contig_" c "\n" $1}' totaluniqseq > uniq.fasta

#### 11. Extract the forward reads
    sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta

#### 12. cdhit assembly. takes some time.
    module load cd-hit/4.6
    
the -c 0.8 in the line below means that you're aligning sequences that have 0.8 similarity and above (sequence similarity threshold - the default is 0.9, so increase it from 0.8 to 0.9)
#do 90, 93, and 95 and see how many contigs we have - run each of these and then record number of contigs
#in the lines below, cdhit is using uniq.F.fasta to create/write (or overwrite) the file called reference.fasta.original. First line listed below takes about 5 mins.

    nohup cd-hit-est -i uniq.F.fasta -o reference.fasta.original -M 0 -T 0 -c 0.8 &>/dev/null &
    grep "^>" reference.fasta.original -c

206818 contigs (grepped below)

    nohup cd-hit-est -i uniq.F.fasta -o reference.fasta.original -M 0 -T 0 -c 0.9 &>/dev/null &
    grep "^>" reference.fasta.original -c
    
349768 contigs

MOVING FORWARD WITH 0.9 MATCH THRESHOLD

### GET RID OF .UNIQ.SEQS FILES!!!

    /working/lgalland/pines_combined

    mkdir bwa
    mv /working/lgalland/pines_combined/reference.fasta.original bwa/
    rm reference.fasta.original.clstr
    
    /working/lgalland/pines_combined/bwa/

Artificial reference is: reference.fasta.original
    
    grep "^>" reference.fasta.original -c
349768 contigs for -c 0.9

#### renaming id lines from contig to scaffold:
    
    sed "s/Contig/scaffold/" reference.fasta.original > pine_ref.fasta
pine_ref.fasta is the reference to use below

    rm reference.fasta.original
    grep ">" pine_ref.fasta -c
349768 contigs in the reference

#### index reference genome:
about 1 minute. makes other ref files (.amb, etc.)
    
    module load bwa/0.7.8
    nohup bwa index -p pine_ref -a is pine_ref.fasta &>/dev/null &
    cd /working/lgalland/pines_combined/
    mv *.gz /working/lgalland/pines_combined/bwa/
    cd /working/lgalland/pines_combined/bwa
    
Unzipping fastqs. Originally did about 5 per minute:
    
    nohup gunzip *fastq.gz &>/dev/null &

## ALREADY REMOVED LOW COVERAGE INDIVIDUALS

    /working/lgalland/pines_combined/bwa
    
#### Calculate the mean number of reads per individual, from your individual fastqs

    grep "^@" -c *.fastq > meanReads_perInd.txt &
    
    R
    meanReads <- read.delim("meanReads_perInd.txt", header=FALSE, sep=":")
    mean(meanReads[,2])
2016325 (this is the mean!)

    quit()

####################################################################################
## 5. reference based assembly using bwa/0.7.5a
####################################################################################

map sequneces to the reference genome using bwa via a perl wrapper to iterate over individuals = multiple fastq files. This perl script runs bwa aln and bwa samse (single end reads). This will produce alignments of the reads for each individual, mapped onto the artificial reference genome.

Running BWA. Note that the script must be modified with the correct reference name (which must be changed in TWO places in the script). Example: the runbwa.pl script references sheep_ref in two places. Change to accurate species. NOTE: We do not have this "sheep_ref" as a file. Rather, it represents the conglomeration of the sheep_ref.amb, .ann, .bwt, .fasta, .pac, and .sa files

Edit distance of 4. This step takes several hours (e.g., 6 hours for 279 pines).
**Working with old bwa 7.5 and samtools 1.3 and bcftools 1.3; following monia and alex's notes. 

    cp /working/lgalland/perl_scripts/runbwa.pl /working/lgalland/pines_combined/bwa/
    nano runbwa.pl
    
    module load bwa/0.7.5a
    nohup bwa index -p pine_ref -a is pine_ref.fasta &>/dev/null &
    nohup perl runbwa.pl *fastq &>/dev/null &

#### convert sam to bam
Takes some time. About 10 minutes per 70 individuals.

    /working/lgalland/pines_combined/bwa/sam_sai
    module load samtools/1.3
    nohup perl /working/lgalland/perl_scripts/sam2bamV1.3.pl *.sam &>/dev/null &

Few minutes for the following:

    module load samtools/1.3
    nohup perl /working/lgalland/perl_scripts/count_assembled.pl *.sorted.bam &>/dev/null &
    
    rm -rf *.sam
    rm -rf *.sai

#### make a plot in R of assembled vs. raw read counts per individual
    scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/pines_combined/bwa/sam_sai/assembled_perind.txt /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/bwa
    
###### In R:
    setwd("/Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/bwa")

	pdf(file="number_of_reads.pdf", width=6, height=6)

	reads<-read.csv("assembled_perind.txt", header=F)
	hist(reads[,2], breaks=15, col="grey", main="")

	dev.off()

####################################################################################
## 6. calling variants 
####################################################################################

    module load bcftools/1.3
    module load samtools/1.3
    
    mv /working/lgalland/pines_combined/bwa/pine_ref* /working/lgalland/pines_combined/bwa/sam_sai/
    /working/lgalland/pines_combined/bwa/sam_sai/

This takes several (2-4) hours. Run as one large chunk of code:

    samtools mpileup -P ILLUMINA --BCF --max-depth 100 --adjust-MQ 50 --min-BQ 20 --min-MQ 20 --skip-indels --output-tags DP,AD --fasta-ref pine_ref.fasta aln*sorted.bam | \
    bcftools call -m --variants-only --format-fields GQ --skip-variants indels | \
    bcftools filter --set-GTs . -i 'QUAL > 19 && FMT/GQ >9' | \
    bcftools view -m 2 -M 2 -v snps --apply-filter "PASS" --output-type v --output-file variants_rawfiltered_9DEC2020.vcf &









# DONE TO HERE. LG 12/8 4:40pm











#### making id file for reheadering:
    sed -s "s/aln_//" assembled_perind.txt | sed -s "s/.sorted.bam//" > pine_ids_col.txt

    module load bcftools/1.3
    module load samtools/1.3

This step is reheadering variants_rawfiltered_19MAR2020.vcf:
    
    nohup bcftools reheader -s pine_ids_col.txt variants_rawfiltered_19MAR2020.vcf -o rehead_variants_rawfiltered_19MAR2020.vcf &>/dev/null & 

####################################################################################
## 7. filtering
####################################################################################

#### common MAF 0.05
VCFtools (0.1.14)

Adam Auton and Anthony Marcketta 2009

Filtering the entire vcf files, at least 90% of individuals have at least one read, and filtering on maf (filtering on loci)
	
	module load vcftools/0.1.14
	module load bcftools/1.3

	vcftools --vcf rehead_variants_rawfiltered_19MAR2020.vcf --out variants_miss10_common --remove-filtered-all --maf 0.05 --max-		missing 0.9 --recode --thin 90 &
** After filtering, kept 279 out of 279 individuals
** After filtering, kept 1584 out of a possible 1743422 Sites 

	vcftools --vcf rehead_variants_rawfiltered_19MAR2020.vcf --out variants_miss20_common --remove-filtered-all --maf 0.05 --max-		missing 0.8 --recode --thin 90 &
** After filtering, kept 279 out of 279 Individuals
** After filtering, kept 5553 out of a possible 1743422 Sites 

	vcftools --vcf rehead_variants_rawfiltered_19MAR2020.vcf --out variants_miss30_common --remove-filtered-all --maf 0.05 --max-		missing 0.7 --recode --thin 90 &
** After filtering, kept 279 out of 279 Individuals
** After filtering, kept 5553 out of a possible 1743422 Sites

	vcftools --vcf rehead_variants_rawfiltered_19MAR2020.vcf --out variants_miss40_common --remove-filtered-all --maf 0.05 --max-		missing 0.6 --recode --thin 90 &
** After filtering, kept 279 out of 279 Individuals
** After filtering, kept 5553 out of a possible 1743422 Sites

#### generate depth of coverage files

	vcftools --vcf variants_miss10_common.recode.vcf --depth -c > variants_miss10_common.txt &
	vcftools --vcf variants_miss20_common.recode.vcf --depth -c > variants_miss20_common.txt &
	vcftools --vcf variants_miss30_common.recode.vcf --depth -c > variants_miss30_common.txt &
	vcftools --vcf variants_miss40_common.recode.vcf --depth -c > variants_miss40_common.txt &

#### generate .mpgl and pntest for each (pntest is genotype likelihood file - use for PCA)

	perl /working/lgalland/perl_scripts/vcf2mpglV1.3TLP.pl variants_miss10_common.recode.vcf
	perl /working/lgalland/perl_scripts/vcf2mpglV1.3TLP.pl variants_miss20_common.recode.vcf
	perl /working/lgalland/perl_scripts/vcf2mpglV1.3TLP.pl variants_miss30_common.recode.vcf
	perl /working/lgalland/perl_scripts/vcf2mpglV1.3TLP.pl variants_miss40_common.recode.vcf

The following writes variants_miss40_common.recode.mpgl, variants_miss30_common.recode.mpgl, etc.

	perl /working/lgalland/perl_scripts/gl2genestV1.3.pl variants_miss10_common.recode.mpgl mean
	perl /working/lgalland/perl_scripts/gl2genestV1.3.pl variants_miss20_common.recode.mpgl mean
	perl /working/lgalland/perl_scripts/gl2genestV1.3.pl variants_miss30_common.recode.mpgl mean
	perl /working/lgalland/perl_scripts/gl2genestV1.3.pl variants_miss40_common.recode.mpgl mean

	module load bcftools/1.3
	module load samtools/1.3

The following cuts the first column and puts it into a new file, which is labeled "good head," but it really isn't. 
	
	cut -d "," -f 1 pine_ids_col.txt > pine_ids_good_head.txt
	
Need to nano the file and add "ind" to the top of the column so it's in proper format:
	
	nano pine_ids_good_head.txt

#### make coverage file
Give it your  mpgl file, your ids file a header (hence "good_head"), and the location of your bam files - RUNNING FOR MISS10, MISS20, MISS30, AND MISS40. (We typically choose one with which we will move forward (e.g., only progressing with miss50), but for radiata, we wanted to see all of the PCAs) (Each takes 5-10 minutes, but is based  on number of individuals. The output file is variants_miss10_common.recode.mpgl_coverage.csv, for example. You can go in another terminal window and verify this file is increasing in size.) Here, for combined pines, we will move forward with miss30.
	
	perl /working/lgalland/perl_scripts/coverage_calc_bam.pl variants_miss30_common.recode.mpgl pine_ids_good_head.txt 			/working/lgalland/pines_combined/bwa/sam_sai/ &
	
	
	
	
	
	
	
	




## variant calling with bcftools 1.9
## bams is a text file with all of the sorted bam files listed, one per line
bcftools mpileup -C 50 -d 250 -f genome.fasta -q 30 -Q 20 -I -b bams -O b -o tpod.bcf

## sometimes I use the -c option, sometimes not, I have mixed feelings, probably would use it
## for podura
bcftools call -v -c -f GQ -p 0.01 -P 0.001 -O v -o tpod.vcf tpod.bcf
