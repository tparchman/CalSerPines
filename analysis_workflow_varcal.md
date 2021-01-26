example of inserting picutre, pdf, etc:
![PCA](https://github.com/tparchman/CalSerPines/blob/master/preliminary_figures/pinesCombined_firstLook_miss30_bySpecies.pdf)

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

Being executed on ponderosa using tapioca pipeline. Commands in bash script, executed as below (11/17/20). This was for two NovaSeq lanes generated in November.

    $ module load fqutils/0.4.1
    $ module load bowtie2/2.2.5
    
    $ bash cleaning_bash.sh &

A second pair of NovaSeq lanes from UTGSAF was processed starting 12/26/20, as below

    $ module load fqutils/0.4.1
    $ module load bowtie2/2.2.5
    
    $ bash cleaning_bash_12_20.sh &

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
    
Parsing BHS_PIMU library :

 	$ perl parse_barcodes768.pl barcodeKey_lib6_bighorns_pines.csv BHS_PIMU.clean.fastq A00 &

Parsing TICR_PIMU library:

	$ perl parse_barcodes768.pl barcodeKey_lib4_timema_pines.csv TICR_PIMU.clean.fastq A00 &

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

    $ cut -f 3 -d "," barcodeKey_lib4_timema_pines.csv | grep "_" > TICR_PIMU_ids_noheader.txt

    $ cut -f 3 -d "," barcodeKey_lib6_bighorns_pines.csv | grep "_" > BHS_PIMU_ids_noheader.txt

Split fastqs by individual

    $ perl splitFastq_universal_regex.pl L1_ids_noheader.txt parsed_S1_11_20.clean.fastq &

    $ perl splitFastq_universal_regex.pl L2_ids_noheader.txt parsed_S2_11_20.clean.fastq &

    $ perl splitFastq_universal_regex.pl TICR_PIMU_ids_noheader.txt parsed_TICR_PIMU.clean.fastq &

    $ perl splitFastq_universal_regex.pl BHS_PIMU_ids_noheader.txt parsed_BHS_PIMU.clean.fastq & 

Zip the parsed*fastq files for now, but delete once patterns and qc are verified:

    $ gzip parsed_S1_11_20.clean.fastq
    $ gzip parsed_S2_11_20.clean.fastq

Total reads for muricata, radiata, and attenuata

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

## STOP!!!

All individuals fastqs need to be backed up on a lab hard drive before continuing. THIS IS VERY IMPORTANT

## Lanie started here: all species
 
    $ cd /working/lgalland/pines_combined/

####################################################################################
## 4. de novo assembly
####################################################################################

### removed individuals with <35M (zipped) of sequence data
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
    rm -rf PM_DM_0241.fastq
    rm -rf PR_MP_0823.fastq

removed 19 individuals. new total = 570 individuals

### Combined pines, de novo assembly

#### 1. gzip files (takes time)
    
    $ nohup gzip *fastq &>/dev/null &

#### 2. make list of fastqs
    
    $ ls *.fastq.gz > namelist

#### 3. one-liner to remove '.fastq.gz' from namelist (instant)
    
    $ sed -i'' -e 's/.fastq.gz//g' namelist
    
#### 4. set variables (run as one chunk command, instant)

    $  AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
    AWK2='!/>/'
    AWK3='!/NNN/'
    PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
    
#### 5. create list of sequences and counts for each individual (few minutes)
    
    $ cat namelist | parallel --no-notice -j 8 "zcat {}.fastq | mawk '$AWK1' | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs" &
    
#### 6. combine files for all individuals (about 15 seconds)
    
    $ cat *.uniq.seqs > uniq.seqs

#### 7. select those sequences that have 4 reads (or however you want to choose - pretty instant, about 1 minute)
    
    $ parallel --no-notice -j 8 mawk -v x=4 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv
    
    $ wc -l uniqCperindv
        16759980 uniq sequences per individual

#### 8. restrict the data by how many individuals have that read
    `NOTE`: the bash below loop just makes a list from 2 to 10

    $ for ((i = 2; i <= 10; i++));
    do
    echo $i >> ufile
    done

    $ cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data
    
    $ more uniqseq.peri.data

        2	5937752
        3	3674369
        4	2709744
        5	2167809
        6	1810265
        7	1550095
        8	1353227
        9	1196393
        10	1068049

    $ rm -rf ufile
    
#### 9. restrict sequences to those that are found in at least 4 inds (few seconds per ind)

    $ mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.4.seqs
    $ wc -l uniq.k.4.c.4.seqs
        2709744 uniq.k.4.c.4.seqs

#### 10. Convert these sequences to fasta format (both pretty instant)
 
    $ cut -f2 uniq.k.4.c.4.seqs > totaluniqseq
    $ mawk '{c= c + 1; print ">Contig_" c "\n" $1}' totaluniqseq > uniq.fasta

#### 11. Extract the forward reads
    
    $ sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta

#### 12. cdhit assembly. takes some time.
   
    $ module load cd-hit/4.6
    
`NOTE`: the -c 0.8 in the line below means that you're aligning sequences that have 0.8 similarity and above (sequence similarity threshold - the default is 0.9, so increase it from 0.8 to 0.9). Do 90, 93, and 95 and see how many contigs we have - run each of these and then record number of contigs.

`NOTE`: in the lines below, `cdhit` is using uniq.F.fasta to create/write (or overwrite) the file called reference.fasta.original. First line listed below takes about 5 mins.

    $ nohup cd-hit-est -i uniq.F.fasta -o reference.fasta.original -M 0 -T 0 -c 0.9 &>/dev/null &

    $ grep "^>" reference.fasta.original -c
        615977 contigs

MOVING FORWARD WITH 0.9 MATCH THRESHOLD

### GET RID OF .UNIQ.SEQS FILES!!!

    $ cd /working/lgalland/pines_combined

    $ mkdir bwa
    $ mv /working/lgalland/pines_combined/reference.fasta.original bwa/
    $ rm reference.fasta.original.clstr
    
    $ cd /working/lgalland/pines_combined/bwa/

Artificial reference is: reference.fasta.original
    
    $ grep "^>" reference.fasta.original -c
        615977 contigs for -c 0.9

Rename id lines from contig to scaffold:
    
    $ sed "s/Contig/scaffold/" reference.fasta.original >pine_ref.fasta

`NOTE`: pine_ref.fasta is the reference to use below

    $ rm reference.fasta.original
    $ grep ">" pine_ref.fasta -c
        615977 contigs in the reference

Index reference genome (about 1 minute. makes other ref files (.amb, etc.))
    
    $ module load bwa/0.7.8
    $ nohup bwa index -p PA_ref -a is pine_ref.fasta &>/dev/null &

    $ cd /working/lgalland/pines_combined/
    $ mv *.gz /working/lgalland/pines_combined/bwa/
    $ cd /working/lgalland/pines_combined/bwa
    
Unzip fastqs. Originally did about 5 per minute:
    
    $ nohup gunzip *fastq.gz &>/dev/null &

    $ cd /working/lgalland/pines_combined/bwa

##########################################
##########################################
##########################################
 # THE FOLLOWING IS NOT IN JAHNER'S FILE YET!!!   
#### Calculate the mean number of reads per individual, from your individual fastqs

    grep "^@" -c *.fastq > meanReads_perInd.txt &
    
    R
    meanReads <- read.delim("meanReads_perInd.txt", header=FALSE, sep=":")
    mean(meanReads[,2])
2590771 (this is the mean!)

    quit()
##########################################
##########################################
##########################################

####################################################################################
## 5. reference based assembly using bwa/0.7.5a
####################################################################################

`NOTE`:  Map sequneces to the reference genome using `bwa` via a perl wrapper to iterate over individuals = multiple fastq files. This perl script runs bwa aln and bwa samse (single end reads). This will produce alignments of the reads for each individual, mapped onto the artificial reference genome.

Running `bwa`. Note that the script must be modified with the correct reference name (which must be changed in TWO places in the script). Example: the `runbwa.pl` script references `sheep_ref` in two places. Change to accurate species. `NOTE`: We do not have this "sheep_ref" as a file. Rather, it represents the conglomeration of the sheep_ref.amb, .ann, .bwt, .fasta, .pac, and .sa files. Edit distance of 4. This step takes several hours (e.g., 6 hours for 279 pines).
**Working with old bwa 7.5 and samtools 1.3 and bcftools 1.3; following monia and alex's notes. 

`NOTE`: `runbwa.pl` needs to be modified for every project (index name, output directory, number of cpus). 

    $ cp /working/lgalland/perl_scripts/runbwa.pl /working/lgalland/pines_combined/bwa/
    $ nano runbwa.pl

 Make index for `bwa` 

    $ module load bwa/0.7.5a
    $ nohup bwa index -p pine_ref -a is pine_ref.fasta &>/dev/null &
 
 Perl wrapper to run `bwa`.    
    
    $ nohup perl runbwa.pl *fastq &>/dev/null &

Convert `sam` to `bam` files. `NOTE`: change number of threads based on current server usage. (Takes some time. About 10 minutes per 70 individuals.)

    $ cd /working/lgalland/pines_combined/bwa/sam_sai
    $ module load samtools/1.3
    $ nohup perl /working/lgalland/perl_scripts/sam2bamV1.3.pl *.sam &>/dev/null &

Clean up a bit

    $ rm -rf *.sam
    $ rm -rf *.sai

Count assembled

    $ module load samtools/1.3
    $ nohup perl /working/lgalland/perl_scripts/count_assembled.pl *.sorted.bam &>/dev/null &
    
Make a plot in R of assembled vs. raw read counts per individual


    $ scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/pines_combined/bwa/sam_sai/assembled_perind.txt /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/bwa
    
`NOTE`: completed in `R`:

    setwd("/Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/bwa")

	pdf(file="number_of_reads.pdf", width=6, height=6)

	reads<-read.csv("assembled_perind.txt", header=F)
	hist(reads[,2], breaks=15, col="grey", main="")

	dev.off()

####################################################################################
## 6. calling variants 
####################################################################################

    $ module load bcftools/1.3
    $ module load samtools/1.3
    
    $ mv /working/lgalland/pines_combined/bwa/pine_ref* /working/lgalland/pines_combined/bwa/sam_sai/
    
    $ cd /working/lgalland/pines_combined/bwa/sam_sai/

The following takes several (2-4) hours. `NOTE`: Run the following lines as one large chunk of code. Be sure to change reference name and output file.

    $ samtools mpileup -P ILLUMINA --BCF --max-depth 100 --adjust-MQ 50 --min-BQ 20 --min-MQ 20 --skip-indels --output-tags DP,AD --fasta-ref pine_ref.fasta aln*sorted.bam | \
    bcftools call -m --variants-only --format-fields GQ --skip-variants indels | \
    bcftools filter --set-GTs . -i 'QUAL > 19 && FMT/GQ >9' | \
    bcftools view -m 2 -M 2 -v snps --apply-filter "PASS" --output-type v --output-file variants_rawfiltered_12JAN2021.vcf &

Number of loci in the vcf

    $ grep -c "^scaffold" variants_rawfiltered_12JAN2021.vcf
        2461844

Make id file for reheadering

    $ sed -s "s/aln_//" assembled_perind.txt | sed -s "s/.sorted.bam//" > pine_ids_col.txt

Reheader variants_rawfiltered_dateHERE.vcf:
    
    $ module load bcftools/1.3
    $ module load samtools/1.3
	$ module load vcftools/0.1.14
    
    $ nohup bcftools reheader -s pine_ids_col.txt variants_rawfiltered_12JAN2021.vcf -o rehead_variants_rawfiltered_12JAN2021.vcf &>/dev/null & 

####################################################################################
## 7. filtering
####################################################################################

Initial round of filtering (just getting a feeling of how filtering parameters might shape the dataset). 570 individuals.

    $ vcftools --vcf rehead_variants_rawfiltered_12JAN2021.vcf --out variants_miss10_maf05 --remove-filtered-all --maf 0.05 --max-missing 0.9 --recode --thin 90
            After filtering, kept 3071 out of a possible 2461844 Sites
    $ vcftools --vcf rehead_variants_rawfiltered_12JAN2021.vcf --out variants_miss20_maf05 --remove-filtered-all --maf 0.05 --max-missing 0.8 --recode --thin 90
	        After filtering, kept 9703 out of a possible 2461844 Sites
    $ vcftools --vcf rehead_variants_rawfiltered_12JAN2021.vcf --out variants_miss30_maf05 --remove-filtered-all --maf 0.05 --max-missing 0.7 --recode --thin 90
	        After filtering, kept 14627 out of a possible 2461844 Sites
    $ vcftools --vcf rehead_variants_rawfiltered_12JAN2021.vcf --out variants_miss40_maf05 --remove-filtered-all --maf 0.05 --max-missing 0.6 --recode --thin 90
	        After filtering, kept 19564 out of a possible 2461844 Sites

`NOTE`: moving forward with `miss30_maf05`

Filter out bad individuals, then refilter

	$ vcftools --vcf variants_miss30_maf05.recode.vcf --missing-indv

`OPTIONAL`: look at distribution of missing data in `R`

    R
    j <- read.delim("out.imiss", header=T)
        dim(j)
        head(j)
    hist(j[,5], breaks=100, col="gray")

`NOTE`: will remove individuals with >50% missing data (N = 17). (As this decimal number decreases, the number of individuals to remove increases. That is, the lower the fraction, the more stringent you are filtering.)

    $ mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
    $ vcftools --vcf rehead_variants_rawfiltered_12JAN2021.vcf --remove lowDP.indv --recode --recode-INFO-all --out rehead_variants_rawfiltered_12JAN2021_noBadInds
        After filtering, kept 553 out of 570 Individuals
        After filtering, kept 2461844 out of a possible 2461844 Sites

Clean up the massive mess of files created above.

	$ rm variants_miss*
    $ gzip variants_rawfiltered_12JAN2021.vcf &
    $ gzip rehead_variants_rawfiltered_12JAN2021.vcf &

Filter the entire raw, NEW vcf file (the one without the bad Inds), at least 70% of individuals have at least one read, and filtering on maf .05 (filtering on loci)

	$ module load vcftools/0.1.14
	$ module load bcftools/1.3

	$ vcftools --vcf rehead_variants_rawfiltered_12JAN2021_noBadInds.recode.vcf --out variants_miss30_maf05 --remove-filtered-all --maf 0.05 --max-missing 0.7 --recode --thin 90
		After filtering, kept 15404 out of a possible 2461844 Sites

Rekill bad inds (missing >50%). 0 inds removed

    $ vcftools --vcf variants_miss30_maf05.recode.vcf --missing-indv
    $ mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
    $ vcftools --vcf variants_miss30_maf05.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out variants_miss30_maf05_noBadInds
			After filtering, kept 553 out of 553 Individuals
            After filtering, kept 15404 out of a possible 15404 Sites

Generate mpgl
	
	$ perl /working/lgalland/perl_scripts/vcf2mpglV1.3TLP.pl variants_miss30_maf05_noBadInds.recode.vcf

Generate pntest genotype likelihood file (for PCA)
	
	$ perl /working/lgalland/perl_scripts/gl2genestV1.3.pl variants_miss30_maf05_noBadInds.recode.mpgl mean

Make a file of IDs. Start with updated vcf file. Cut the first column and puts it into a new file, which is labeled "good head."

    $ vcftools --vcf variants_miss30_maf05_noBadInds.recode.vcf --missing-indv 
    $ cut -f 1 out.imiss > pine_ids_553.txt
    $ sed "s/INDV/ind/" pine_ids_553.txt | sed "s/aln_//g" | sed "s/.sorted.bam//g" > pine_ids_553_good_head.txt

Make other ID and pop formats for PCA. 
	
	$ cp pine_ids_553_good_head.txt pine_ids_553_good_noHead.txt
	$ nano pine_ids_553_good_noHead.txt
		remove the "ind" at the top for proper PCA format. The following tells awk that the delimiter is a comma, and to grab and print the first column, and then print it into pine_ids_col.txt (instead of just printing to screen)

The following tells awk that the delimiter is a comma, and to grab and print the first column, and then print it into pine_ids_col.txt (instead of just printing to screen)

	$ cp pine_ids_553_good_noHead.txt pine_ids_553_good_noHead_original.txt
	$ awk -F "," '{print $1}' pine_ids_553_good_noHead.txt > pine_553_ids_col.txt

Then, `sed` out the parts you don't want for your pops file (removing the "first term" in your ID file, in this case). (The following worked here, which was in format XX(species)_XX(population)_XXXX(ind).) 

	$ sed -s "s/[A-Z][A-Z]_//" pine_553_ids_col.txt > pine_pops3.txt
	$ sed -s "s/_[0-9]*//" pine_pops3.txt > pine_553_pops.txt

***These are the two files you need for initial verification PCA (species_ids_col.txt and species_pops.txt). file `species_num_ids_col.txt` is a single column format, no header, with individuals listed as HT_MR_0692, etc. File `species_num_pops.txt` is a single column format, no header, with lines (in order) as CN, CN, CN, CN, CN, AB, AB, AB, AB.	

`scp` files to laptop for initial look at PCAs

	$ scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/pines_combined/bwa/sam_sai/pntest_mean_variants_miss30_maf05_noBadInds.recode.txt /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/bwa/PCA 

	$ scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/pines_combined/bwa/sam_sai/pine_553_ids_col.txt /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/bwa/PCA 

	$ scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/pines_combined/bwa/sam_sai/pine_553_pops.txt /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/bwa/PCA 

Initial look at PCAs to make sure data is okay, done in R.  

`NOTE`: See new_gl_attempt2.R in /Users/lanie/lanie/PhD/genomics/pines/attenuata/bwa/PCA for the highly edited PCAs where populations are plotted in helpful order, and individuals are removed.

	setwd("/Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/bwa/PCA")

	#miss 30 
	read.table("pntest_mean_variants_miss30_maf05_noBadInds.recode.txt", header=F)->gl
	read.table("pine_553_ids_col.txt", header=F)->ids
	read.table("pine_553_pops.txt", header=F)->pops
	
    # If you need to remove individuals after seeing the PCA, do the following to remove rows. Just don't forget to change the number of individuals later in the script.
        ids <- ids[-c(354, 143:144, 480), ]
        pops <- pops[-c(354, 143:144, 480), ]
        gl <- gl[ , -c(354, 143:144, 480)]

        write.table(ids, file="ids_new.txt", sep=" ", row.names=F, col.names=F , quote=F)
        write.table(pops, file="pops_new.txt", sep=" ", row.names=F, col.names=F , quote=F)
        write.table(gl, file="gl_new.txt", sep=" ", row.names=F, col.names=F , quote=F)

        read.table("gl_new.txt", header=F) -> gl
        read.table("ids_new.txt", header=F) -> ids
        read.table("pops_new.txt", header=F) -> pops


    t(gl)->tgl
	cbind(ids, pops, tgl)->tidsgl
	write.table(tidsgl, file="new_pine_gl_matrix_miss30_maf05.txt", sep=" ", row.names=F, col.names=F , quote=F)

	# Now, the files are ready for PCA
	# miss 30

	miss30 <- read.delim("new_pine_gl_matrix_miss30_maf05.txt", header=FALSE, sep=" ")
	miss30[1:10,1:10]
	dim(miss30) #549 * 15406 ####two more than loci because first column is whole ID identifier, second column is population identifier

	g30 <- t(miss30[,3:15406])
	dim(g30) # 15404 * 549
	g30[1:10,1:10]

	gmn30 <- apply(g30, 1, mean, na.rm=TRUE)
	gmnmat30 <- matrix(gmn30, nrow=15404, ncol=549)
	gprime30 <- g30 - gmnmat30
	gcovarmat30 <- matrix(NA, nrow=549, ncol=549)

	for (i in 1:549)
	{
  		for (j in 1:549)
  		{
  		if (i==j)
  		{
  		gcovarmat30[i,j] <- cov(gprime30[,i], gprime30[,j], use="pairwise.complete.obs")
  		}
  		else
  		{
  		gcovarmat30[i,j] <- cov(gprime30[,i], gprime30[,j], use="pairwise.complete.obs")
  		gcovarmat30[j,i] <- gcovarmat30[i,j]
  		}
  		}
	}

	pcgcov30<-prcomp(x=gcovarmat30,center=TRUE,scale=FALSE)
	imp30 <- summary(pcgcov30)
	summary(pcgcov30)

    pcgcov30$x[,1]
    
	# Plotting PCs, miss30

	colors <-c(
  		"#b34469",
        "#4fb848",
        "#9557cf",
        "#8fb63d",
        "#5a72e0",
        "#c2b43a",
        "#cf5dc8",
        "#52c47c",
        "#d84491",
        "#5b872a",
        "#6e51a3",
        "#dc9130",
        "#5a74b9",
        "#e46332",
        "#3abec8",
        "#d13736",
        "#53c2a3",
        "#db3963",
        "#69ab6e",
        "#9d5096",
        "#35773e",
        "#b68edf",
        "#947c2b",
        "#60a4da",
        "#ae5020",
        "#308266",
        "#bc4c49",
        "#9fab65",
        "#dc87b7",
        "#5f6726",
        "#924a67",
        "#d3a667",
        "#d97b7e",
        "#965e30",
        "#e48e69")

	# PC 1v2, miss30, independent colors

	# the plot function in the line below sets the axes right off the bat
	# so, the pcgcov30$x[,1] refers to PC 1, and pcgcov30$x[,2] refers to PC 2

	plot(pcgcov30$x[,1], pcgcov30$x[,2], type="n", main = "Pines combined, first look miss30", xlab=paste("PC1 (",(imp30$importance[,1][[2]]*100), "% )", sep=""), ylab=paste("PC2 (",(imp30$importance[,2][[2]]*100), "% )", sep=""), cex.lab=1.2)

	# the blue 1 and 2 in the points function refer to the PCs
	# must change these to 2 and 3 for PCs 2 v 3, for example!

    # P. muricata
    points(pcgcov30$x[which(pops=="DC"),1], pcgcov30$x[which(pops=="DC"), 2], pch=21, bg=colors[1], cex=1)
    points(pcgcov30$x[which(pops=="FR"),1], pcgcov30$x[which(pops=="FR"), 2], pch=21, bg=colors[2], cex=1)
    points(pcgcov30$x[which(pops=="NR"),1], pcgcov30$x[which(pops=="NR"), 2], pch=21, bg=colors[3], cex=1)
    points(pcgcov30$x[which(pops=="PA"),1], pcgcov30$x[which(pops=="PA"), 2], pch=21, bg=colors[4], cex=1)
    points(pcgcov30$x[which(pops=="PP"),1], pcgcov30$x[which(pops=="PP"), 2], pch=21, bg=colors[5], cex=1)
    points(pcgcov30$x[which(pops=="SP"),1], pcgcov30$x[which(pops=="SP"), 2], pch=21, bg=colors[6], cex=1)
    points(pcgcov30$x[which(pops=="SR"),1], pcgcov30$x[which(pops=="SR"), 2], pch=21, bg=colors[7], cex=1)
    points(pcgcov30$x[which(pops=="PR"),1], pcgcov30$x[which(pops=="PR"), 2], pch=21, bg=colors[8], cex=1)
    points(pcgcov30$x[which(pops=="CH"),1], pcgcov30$x[which(pops=="CH"), 2], pch=21, bg=colors[9], cex=1)
    points(pcgcov30$x[which(pops=="CP"),1], pcgcov30$x[which(pops=="CP"), 2], pch=21, bg=colors[10], cex=1)
    points(pcgcov30$x[which(pops=="DM"),1], pcgcov30$x[which(pops=="DM"), 2], pch=21, bg=colors[11], cex=1)
    points(pcgcov30$x[which(pops=="LO"),1], pcgcov30$x[which(pops=="LO"), 2], pch=21, bg=colors[12], cex=1)
    points(pcgcov30$x[which(pops=="PB"),1], pcgcov30$x[which(pops=="PB"), 2], pch=21, bg=colors[13], cex=1)
    points(pcgcov30$x[which(pops=="RR"),1], pcgcov30$x[which(pops=="RR"), 2], pch=21, bg=colors[14], cex=1)

    # P. attenuata
	points(pcgcov30$x[which(pops=="AH"),1], pcgcov30$x[which(pops=="AH"), 2], pch=21, bg=colors[15], cex=1)
	points(pcgcov30$x[which(pops=="AL"),1], pcgcov30$x[which(pops=="AL"), 2], pch=21, bg=colors[16], cex=1)
	points(pcgcov30$x[which(pops=="BS"),1], pcgcov30$x[which(pops=="BS"), 2], pch=21, bg=colors[17], cex=1)
	points(pcgcov30$x[which(pops=="CG"),1], pcgcov30$x[which(pops=="CG"), 2], pch=21, bg=colors[18], cex=1)
	points(pcgcov30$x[which(pops=="LA"),1], pcgcov30$x[which(pops=="LA"), 2], pch=21, bg=colors[19], cex=1)
	points(pcgcov30$x[which(pops=="LS"),1], pcgcov30$x[which(pops=="LS"), 2], pch=21, bg=colors[20], cex=1)
	points(pcgcov30$x[which(pops=="OC"),1], pcgcov30$x[which(pops=="OC"), 2], pch=21, bg=colors[21], cex=1)
	points(pcgcov30$x[which(pops=="PF"),1], pcgcov30$x[which(pops=="PF"), 2], pch=21, bg=colors[22], cex=1)
	points(pcgcov30$x[which(pops=="SB"),1], pcgcov30$x[which(pops=="SB"), 2], pch=21, bg=colors[23], cex=1)
	points(pcgcov30$x[which(pops=="SH"),1], pcgcov30$x[which(pops=="SH"), 2], pch=21, bg=colors[24], cex=1)
	points(pcgcov30$x[which(pops=="SL"),1], pcgcov30$x[which(pops=="SL"), 2], pch=21, bg=colors[25], cex=1)
	points(pcgcov30$x[which(pops=="ST"),1], pcgcov30$x[which(pops=="ST"), 2], pch=21, bg=colors[26], cex=1)
	points(pcgcov30$x[which(pops=="YA"),1], pcgcov30$x[which(pops=="YA"), 2], pch=21, bg=colors[27], cex=1)
	points(pcgcov30$x[which(pops=="YB"),1], pcgcov30$x[which(pops=="YB"), 2], pch=21, bg=colors[28], cex=1)

    # P. radiata
    points(pcgcov30$x[which(pops=="CM"),1], pcgcov30$x[which(pops=="CM"), 2], pch=21, bg=colors[29], cex=1)
    points(pcgcov30$x[which(pops=="CN"),1], pcgcov30$x[which(pops=="CN"), 2], pch=21, bg=colors[30], cex=1)
    points(pcgcov30$x[which(pops=="CS"),1], pcgcov30$x[which(pops=="CS"), 2], pch=21, bg=colors[31], cex=1)
    points(pcgcov30$x[which(pops=="GU"),1], pcgcov30$x[which(pops=="GU"), 2], pch=21, bg=colors[32], cex=1)
    points(pcgcov30$x[which(pops=="MR"),1], pcgcov30$x[which(pops=="MR"), 2], pch=21, bg=colors[33], cex=1)
    points(pcgcov30$x[which(pops=="SC"),1], pcgcov30$x[which(pops=="SC"), 2], pch=21, bg=colors[34], cex=1)
	points(pcgcov30$x[which(pops=="MP"),1], pcgcov30$x[which(pops=="MP"), 2], pch=21, bg=colors[35], cex=1)

    # Morphological hybrids are listed along their sites. All are part of Santa Cruz (SC) or MR (Monterey) clusters

    legend("bottomright", legend=c("Diablo Canyon PIMU", "Fort Ross PIMU", "Navarro River PIMU", "Point Arena PIMU", "Patrick Point PIMU", "Salt Point PIMU", "Sea Ranch PIMU", "Point Reyes PIMU", "China Pines SCI PIMU", "Christy Pines SCI PIMU", "Monterey Peninsula (DM) PIMU", "Lompoc PIMU", "Pelican Bay SCI PIMU", "Ridge Road PIMU", "Auburn High PIAT", "Auburn Low PIAT", "Big Sur PIAT", "Cuesta Grade(SLO) PIAT", "Los Angeles PIAT", "Lake Shasta PIAT", "Orange County PIAT", "Panther Flat (OR border) PIAT", "San Bernardino PIAT", "Santa Cruz High PIAT", "Santa Cruz Low PIAT", "Santa Cruz Top PIAT", "Yosemite A PIAT", "Yosemite B PIAT", "Cambria PIRA", "Cedros North PIRA", "Cedros South PIRA", "Guadalupe PIRA", "Monterey PIRA", "Santa Cruz PIRA", "Monterey Peninsula PIRA"), pch=c(16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16), ncol=2, col=colors[1:35], cex=.38)

	# save plot window as pinesCombined_firstLook_miss30_maf05.pdf

    ####### PC 1v2, miss30, grouped on species

    colors2 <-c(
  		"#1694e8",
		"#f7b84f",
		"#198c49")

    plot(pcgcov30$x[,1], pcgcov30$x[,2], type="n", main = "Pines combined, FirstLook miss 30, by species", xlab=paste("PC1 (",(imp30$importance[,1][[2]]*100), "% )", sep=""), ylab=paste("PC2 (",(imp30$importance[,2][[2]]*100), "% )", sep=""), cex.lab=1.2)

    # P. muricata
    points(pcgcov30$x[which(pops=="DC"),1], pcgcov30$x[which(pops=="DC"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="FR"),1], pcgcov30$x[which(pops=="FR"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="NR"),1], pcgcov30$x[which(pops=="NR"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="PA"),1], pcgcov30$x[which(pops=="PA"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="PP"),1], pcgcov30$x[which(pops=="PP"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="SP"),1], pcgcov30$x[which(pops=="SP"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="SR"),1], pcgcov30$x[which(pops=="SR"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="PR"),1], pcgcov30$x[which(pops=="PR"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="CH"),1], pcgcov30$x[which(pops=="CH"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="CP"),1], pcgcov30$x[which(pops=="CP"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="DM"),1], pcgcov30$x[which(pops=="DM"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="LO"),1], pcgcov30$x[which(pops=="LO"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="PB"),1], pcgcov30$x[which(pops=="PB"), 2], pch=21, bg=colors2[1], cex=1)
    points(pcgcov30$x[which(pops=="RR"),1], pcgcov30$x[which(pops=="RR"), 2], pch=21, bg=colors2[1], cex=1)

    # P. attenuata:
    points(pcgcov30$x[which(pops=="AH"),1], pcgcov30$x[which(pops=="AH"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="AL"),1], pcgcov30$x[which(pops=="AL"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="BS"),1], pcgcov30$x[which(pops=="BS"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="CG"),1], pcgcov30$x[which(pops=="CG"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="LA"),1], pcgcov30$x[which(pops=="LA"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="LS"),1], pcgcov30$x[which(pops=="LS"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="OC"),1], pcgcov30$x[which(pops=="OC"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="PF"),1], pcgcov30$x[which(pops=="PF"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="SB"),1], pcgcov30$x[which(pops=="SB"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="SH"),1], pcgcov30$x[which(pops=="SH"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="SL"),1], pcgcov30$x[which(pops=="SL"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="ST"),1], pcgcov30$x[which(pops=="ST"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="YA"),1], pcgcov30$x[which(pops=="YA"), 2], pch=21, bg=colors2[3], cex=1)
    points(pcgcov30$x[which(pops=="YB"),1], pcgcov30$x[which(pops=="YB"), 2], pch=21, bg=colors2[3], cex=1)

    # P. radiata
    points(pcgcov30$x[which(pops=="CM"),1], pcgcov30$x[which(pops=="CM"), 2], pch=21, bg=colors2[2], cex=1)
    points(pcgcov30$x[which(pops=="CN"),1], pcgcov30$x[which(pops=="CN"), 2], pch=21, bg=colors2[2], cex=1)
    points(pcgcov30$x[which(pops=="CS"),1], pcgcov30$x[which(pops=="CS"), 2], pch=21, bg=colors2[2], cex=1)
    points(pcgcov30$x[which(pops=="GU"),1], pcgcov30$x[which(pops=="GU"), 2], pch=21, bg=colors2[2], cex=1)
    points(pcgcov30$x[which(pops=="MR"),1], pcgcov30$x[which(pops=="MR"), 2], pch=21, bg=colors2[2], cex=1)
    points(pcgcov30$x[which(pops=="SC"),1], pcgcov30$x[which(pops=="SC"), 2], pch=21, bg=colors2[2], cex=1)
	points(pcgcov30$x[which(pops=="MP"),1], pcgcov30$x[which(pops=="MP"), 2], pch=21, bg=colors2[2], cex=1)

    legend("bottomright", legend=c("P. muricata", "P. radiata", "P. attenuata"), pch=c(16,16,16), ncol=1, col=colors2[1:3], cex=1)

    # save window as pinesCombined_firstLook_bySpecies.pdf
	# copy and paste back to notes file

![PCA_pinesCombined](https://github.com/tparchman/CalSerPines/blob/master/preliminary_figures/pinesCombined_firstLook_miss30_maf05.pdf) 
![PCA_pinesCombined_groupedBySpecies](https://github.com/tparchman/CalSerPines/blob/master/preliminary_figures/pinesCombined_firstLook_miss30_bySpecies.pdf) 

    ############### quick and dirty to plot only muricata:
    
    #miss 30 
    read.table("pntest_mean_variants_miss30_maf05_noBadInds.recode.txt", header=F)->gl
    read.table("pine_553_ids_col.txt", header=F)->ids
    read.table("pine_553_pops.txt", header=F)->pops

    # If you need to remove individuals after seeing the PCA, do the following to remove rows. Just don't forget to change the number of individuals later in the script.
    ids <- ids[-c(1:189, 233, 335, 354, 426:553), ]
    pops <- pops[-c(1:189, 233, 335, 354, 426:553), ]
    gl <- gl[ , -c(1:189, 233, 335, 354, 426:553)]

    write.table(ids, file="ids_new.txt", sep=" ", row.names=F, col.names=F , quote=F)
    write.table(pops, file="pops_new.txt", sep=" ", row.names=F, col.names=F , quote=F)
    write.table(gl, file="gl_new.txt", sep=" ", row.names=F, col.names=F , quote=F)

    read.table("gl_new.txt", header=F) -> gl
    read.table("ids_new.txt", header=F) -> ids
        read.table("pops_new.txt", header=F) -> pops
    dim(gl)
    dim(ids)
    dim(pops)

    t(gl)->tgl
    cbind(ids, pops, tgl)->tidsgl
    write.table(tidsgl, file="pine_gl_matrix_new.txt", sep=" ", row.names=F, col.names=F , quote=F)

    # Now, the files are ready for PCA
    # miss 30

    miss30 <- read.delim("pine_gl_matrix_new.txt", header=FALSE, sep=" ")
    miss30[1:10,1:10]
    dim(miss30) #233 * 15406 ####two more than loci because first column is whole ID identifier, second column  is population identifier

    g30 <- t(miss30[,3:15406])
    dim(g30) # 15404 * 233
    g30[1:10,1:10]

    gmn30 <- apply(g30, 1, mean, na.rm=TRUE)
    gmnmat30 <- matrix(gmn30, nrow=15404, ncol=233)
    gprime30 <- g30 - gmnmat30
    gcovarmat30 <- matrix(NA, nrow=233, ncol=233)

    for (i in 1:233)
        {
        for (j in 1:233)
        {
        if (i==j)
        {
        gcovarmat30[i,j] <- cov(gprime30[,i], gprime30[,j], use="pairwise.complete.obs")
         }
        else
         {
          gcovarmat30[i,j] <- cov(gprime30[,i], gprime30[,j], use="pairwise.complete.obs")
          gcovarmat30[j,i] <- gcovarmat30[i,j]
        }
        }
    }

    pcgcov30<-prcomp(x=gcovarmat30,center=TRUE,scale=FALSE)
    imp30 <- summary(pcgcov30)
    summary(pcgcov30)

    pcgcov30$x[,1]

    # Plotting PCs, miss30

    colors <-c(
        "#b34469",
        "#4fb848",
        "#9557cf",
        "#8fb63d",
        "#5a72e0",
        "#c2b43a",
        "#cf5dc8",
        "#52c47c",
        "#d84491",
        "#5b872a",
        "#6e51a3",
        "#dc9130",
        "#5a74b9",
        "#e46332")

    # PC 1v2, miss30, independent colors

    # the plot function in the line below sets the axes right off the bat
    # so, the pcgcov30$x[,1] refers to PC 1, and pcgcov30$x[,2] refers to PC 2

    plot(pcgcov30$x[,1], pcgcov30$x[,2], type="n", main = "QuickAndDirty muricata, first look miss30", xlab=paste("PC1 (",(imp30$importance[,1][[2]]*100), "% )", sep=""), ylab=paste("PC2 (",(imp30$importance[,2][[2]]*100), "% )", sep=""), cex.lab=1.2)

    # the blue 1 and 2 in the points function refer to the PCs
    # must change these to 2 and 3 for PCs 2 v 3, for example!

    # P. muricata
    points(pcgcov30$x[which(pops=="DC"),1], pcgcov30$x[which(pops=="DC"), 2], pch=21, bg=colors[1], cex=1)
    points(pcgcov30$x[which(pops=="FR"),1], pcgcov30$x[which(pops=="FR"), 2], pch=21, bg=colors[2], cex=1)
    points(pcgcov30$x[which(pops=="NR"),1], pcgcov30$x[which(pops=="NR"), 2], pch=21, bg=colors[3], cex=1)
    points(pcgcov30$x[which(pops=="PA"),1], pcgcov30$x[which(pops=="PA"), 2], pch=21, bg=colors[4], cex=1)
    points(pcgcov30$x[which(pops=="PP"),1], pcgcov30$x[which(pops=="PP"), 2], pch=21, bg=colors[5], cex=1)
    points(pcgcov30$x[which(pops=="SP"),1], pcgcov30$x[which(pops=="SP"), 2], pch=21, bg=colors[6], cex=1)
    points(pcgcov30$x[which(pops=="SR"),1], pcgcov30$x[which(pops=="SR"), 2], pch=21, bg=colors[7], cex=1)
    points(pcgcov30$x[which(pops=="PR"),1], pcgcov30$x[which(pops=="PR"), 2], pch=21, bg=colors[8], cex=1)
    points(pcgcov30$x[which(pops=="CH"),1], pcgcov30$x[which(pops=="CH"), 2], pch=21, bg=colors[9], cex=1)
    points(pcgcov30$x[which(pops=="CP"),1], pcgcov30$x[which(pops=="CP"), 2], pch=21, bg=colors[10], cex=1)
    points(pcgcov30$x[which(pops=="DM"),1], pcgcov30$x[which(pops=="DM"), 2], pch=21, bg=colors[11], cex=1)
    points(pcgcov30$x[which(pops=="LO"),1], pcgcov30$x[which(pops=="LO"), 2], pch=21, bg=colors[12], cex=1)
    points(pcgcov30$x[which(pops=="PB"),1], pcgcov30$x[which(pops=="PB"), 2], pch=21, bg=colors[13], cex=1)
    points(pcgcov30$x[which(pops=="RR"),1], pcgcov30$x[which(pops=="RR"), 2], pch=21, bg=colors[14], cex=1)

    legend("bottomleft", legend=c("Diablo Canyon PIMU", "Fort Ross PIMU", "Navarro River PIMU", "Point Arena PIMU", "Patrick Point PIMU", "Salt Point PIMU", 
    "Sea Ranch PIMU", "Point Reyes PIMU", "China Pines SCI PIMU", "Christy Pines SCI PIMU", "Monterey Peninsula (DM) PIMU", "Lompoc PIMU", "Pelican Bay SCI PIMU", 
    "Ridge Road SCI PIMU"), pch=c(16,16,16,16,16,16,16,16,16,16,16,16,16,16), ncol=2, col=colors[1:14], cex=.8)
                               
    # save plot window as muricata_quickAndDirty_miss30_maf05.pdf

![PCA_muricata_quickAndDirty](https://github.com/tparchman/CalSerPines/blob/master/preliminary_figures/muricata_quickAndDirty_miss30_maf05.pdf) 

    ###### finally, quick and dirty to look at island populations only

    #miss 30 
    read.table("gl_new.txt", header=F)->gl
    read.table("ids_new.txt", header=F)->ids
    read.table("pops_new.txt", header=F)->pops

    # If you need to remove individuals after seeing the PCA, do the following to remove rows. Just don't forget to change the number of individuals later in the script.
    ids <- ids[-c(36, 39:125, 145:176, 196:233), ]
    pops <- pops[-c(36, 39:125, 145:176, 196:233), ]
    gl <- gl[ , -c(36, 39:125, 145:176, 196:233)]

    write.table(ids, file="ids_islandOnly.txt", sep=" ", row.names=F, col.names=F , quote=F)
    write.table(pops, file="pops_islandOnly.txt", sep=" ", row.names=F, col.names=F , quote=F)
    write.table(gl, file="gl_islandOnly.txt", sep=" ", row.names=F, col.names=F , quote=F)

    read.table("gl_islandOnly.txt", header=F) -> gl
    read.table("ids_islandOnly.txt", header=F) -> ids
        read.table("pops_islandOnly.txt", header=F) -> pops
    dim(gl)
    dim(ids)
    dim(pops)

    t(gl)->tgl
    cbind(ids, pops, tgl)->tidsgl
    write.table(tidsgl, file="pine_gl_matrix_islandOnly.txt", sep=" ", row.names=F, col.names=F , quote=F)

    # Now, the files are ready for PCA
    # miss 30

    miss30 <- read.delim("pine_gl_matrix_islandOnly.txt", header=FALSE, sep=" ")
    miss30[1:10,1:10]
    dim(miss30) #75 * 15406 ####two more than loci because first column is whole ID identifier, second column  is population identifier

    g30 <- t(miss30[,3:15406])
    dim(g30) # 15404 * 75
    g30[1:10,1:10]

    gmn30 <- apply(g30, 1, mean, na.rm=TRUE)
    gmnmat30 <- matrix(gmn30, nrow=15404, ncol=75)
    gprime30 <- g30 - gmnmat30
    gcovarmat30 <- matrix(NA, nrow=75, ncol=75)

    for (i in 1:75)
        {
        for (j in 1:75)
        {
        if (i==j)
        {
        gcovarmat30[i,j] <- cov(gprime30[,i], gprime30[,j], use="pairwise.complete.obs")
         }
        else
         {
          gcovarmat30[i,j] <- cov(gprime30[,i], gprime30[,j], use="pairwise.complete.obs")
          gcovarmat30[j,i] <- gcovarmat30[i,j]
        }
        }
    }

    pcgcov30<-prcomp(x=gcovarmat30,center=TRUE,scale=FALSE)
    imp30 <- summary(pcgcov30)
    summary(pcgcov30)

    pcgcov30$x[,1]

    # Plotting PCs, miss30

    colors <-c(
        "#d84491",
        "#5b872a",
        "#5a74b9",
        "#e46332")

    # PC 1v2, miss30, independent colors

    # the plot function in the line below sets the axes right off the bat
    # so, the pcgcov30$x[,1] refers to PC 1, and pcgcov30$x[,2] refers to PC 2

    plot(pcgcov30$x[,1], pcgcov30$x[,2], type="n", main = "QuickAndDirty muricata, SCI only", xlab=paste("PC1 (",(imp30$importance[,1][[2]]*100), "% )", sep=""), ylab=paste("PC2 (",(imp30$importance[,2][[2]]*100), "% )", sep=""), cex.lab=1.2)

    # the blue 1 and 2 in the points function refer to the PCs
    # must change these to 2 and 3 for PCs 2 v 3, for example!

    # P. muricata, Santa Cruz Island only
    points(pcgcov30$x[which(pops=="CH"),1], pcgcov30$x[which(pops=="CH"), 2], pch=21, bg=colors[9], cex=1)
    points(pcgcov30$x[which(pops=="CP"),1], pcgcov30$x[which(pops=="CP"), 2], pch=21, bg=colors[10], cex=1)
    points(pcgcov30$x[which(pops=="PB"),1], pcgcov30$x[which(pops=="PB"), 2], pch=21, bg=colors[13], cex=1)
    points(pcgcov30$x[which(pops=="RR"),1], pcgcov30$x[which(pops=="RR"), 2], pch=21, bg=colors[14], cex=1)

    legend("bottomleft", legend=c("China Pines SCI PIMU", "Christy Pines SCI PIMU", "Pelican Bay SCI PIMU", 
    "Ridge Road SCI PIMU"), pch=c(16,16,16,16), ncol=2, col=colors[1:4], cex=.8)
                               
    # save plot window as muricata_quickAndDirty_islandOnly.pdf

![PCA_muricata_quickAndDirty_islandOnly](https://github.com/tparchman/CalSerPines/blob/master/preliminary_figures/muricata_quickAndDirty_islandOnly.pdf) 









# DONE TO HERE








 
Make coverage file. Give it your  mpgl file, your ids file, and the location of your `bam` files. The output file is *.recode.mpgl_coverage.csv. You can go in another terminal window and verify this file is increasing in size.
	
	perl /working/lgalland/perl_scripts/coverage_calc_bam.pl variants_miss30_maf04_noBadInds.recode.mpgl PA_ids_188_good_head.txt /working/lgalland/attenuata/bwa/sam_sai/ &
	
Need to filter out over-assembled loci. Start by cutting the scaffold names column from the original file into the new file, for use in R.    

    $ cut -d " " -f 1 variants_miss30_maf04_noBadInds.recode.mpgl > loc_names_15495.txt

Move files for use in R.

    $ scp lgalland@134.197.63.151:/working/lgalland/attenuata/bwa/sam_sai/variants_miss30_maf04_noBadInds.recode.mpgl_coverage.csv /Users/lanie/lanie/PhD/genomics/pines/attenuata

    $ scp lgalland@134.197.63.151:/working/lgalland/attenuata/bwa/sam_sai/loc_names_15495.txt /Users/lanie/lanie/PhD/genomics/pines/attenuata

Identify over-assembled loci in R

    R
    
    setwd("/Users/lanie/lanie/PhD/genomics/pines/attenuata")

    dat <- read.csv("variants_miss30_maf04_noBadInds.recode.mpgl_coverage.csv", header=F)
        dim(dat) #188 * 15496 - yes, because first column is ID
        dat[1:10,1:10]

    loc_names <- read.delim("loc_names_15495.txt", header=F)
        head(loc_names)
        dim(loc_names) #15495 * 1

    dat_noname <- dat[,-1] # rip off the first column, which is ID
        dim(dat_noname) #188 * 15495
        dat_noname[1:10,1:10]

    #RERUN THIS FOLLOWING VECTOR EVERY TIME YOU CHANGE PARAMETERS, BECAUSE THE FOR LOOP BELOW APPENDS AND AMMENDS IT
    #start with 10, 12, 15, 20, 25. If no tail is visible, increase to 25, 30, 35, 40, 45

    avg_15495 <- vector()
    in_out_15495_20 <- vector()
    in_out_15495_25 <- vector()
    in_out_15495_30 <- vector()
    in_out_15495_35 <- vector()
    in_out_15495_40 <- vector()

    for (i in 1:15495)
        {
        avg <- mean(dat_noname[,i])
        avg_15495 <- append(avg_15495, avg)
  
        if (avg <= 20)	{ in_out_15495_20 <- append(in_out_15495_20, 1) }
        else			{ in_out_15495_20 <- append(in_out_15495_20, 0) }
  
        if (avg <= 25)	{ in_out_15495_25 <- append(in_out_15495_25, 1) }
        else			{ in_out_15495_25 <- append(in_out_15495_25, 0) }
  
        if (avg <= 30)	{ in_out_15495_30 <- append(in_out_15495_30, 1) }
        else			{ in_out_15495_30 <- append(in_out_15495_30, 0) }
  
        if (avg <= 35)	{ in_out_15495_35 <- append(in_out_15495_35, 1) }
        else			{ in_out_15495_35 <- append(in_out_15495_35, 0) }
  
        if (avg <= 40)	{ in_out_15495_40 <- append(in_out_15495_40, 1) }
        else			{ in_out_15495_40 <- append(in_out_15495_40, 0) }
        }

    hist(avg_15495)
    hist(avg_15495, xlim=c(0,80), breaks=2000) # can change the scale here
    
    # if you look at this histogram and don't see a tail, change the numbers above (10, 12, 15, 20) to other values in the for loop

    ##choose the value at which the numbers STOP declining.

    sum(in_out_15495_20) 
        ## 12612
    sum(in_out_15495_25) 
        ## 13172
    sum(in_out_15495_30) 
        ## 13544
    sum(in_out_15495_35) 	
        ## 13815
    sum(in_out_15495_40)  ## choosing to kill all locs with mean cov/ind >= 40
        ## 14040


    sub_40 <- dat_noname[,in_out_15495_40==1]
        dim(sub_40)
        # 188 x 14040    
    sub_40_avg <- subset(avg_15495, in_out_15495_40==1)	

    kill_locs <- subset(loc_names, in_out_15495_40==0)
        dim(kill_locs)
        # 1455 x 1 (because 15495 - 1455 = 14040, or the number of loci remaining after we remove loci) 
        head(kill_locs)

    write.table(kill_locs, file="high_cov_loc_list_to_be_removed.txt", quote=F, row.names=F, col.names=F)

    #### copy and paste all of this back to the notes file, to keep record of everything that was done!

Time to filter out these over-assembled loci 

    $ scp /Users/lanie/lanie/PhD/genomics/pines/attenuata/high_cov_loc_list_to_be_removed.txt lgalland@134.197.63.151:/working/lgalland/attenuata/bwa/sam_sai

    $ sed "s/:/\t/" high_cov_loc_list_to_be_removed.txt > high_cov_loc_list_to_be_removed_tabdelim.txt

    $ module load vcftools/0.1.14 

    $ vcftools --vcf variants_miss30_maf04_noBadInds.recode.vcf --exclude-positions high_cov_loc_list_to_be_removed_tabdelim.txt --recode --recode-INFO-all --out variants_miss30_maf04_noBadInds_noHighCov
	    ## After filtering, kept 14040 out of a possible 15495 Sites

Generate NEW depth of coverage files, with no bad inds, no over-assembled loci
    $ vcftools --vcf variants_miss30_maf04_noBadInds_noHighCov.recode.vcf --depth -c > variants_miss30_maf04_noBadInds_noHighCov.txt

Create mpgl file NEW .mpgl file (no bad inds, no over-assembled loci)

    $ perl /working/lgalland/perl_scripts/vcf2mpglV1.3TLP.pl variants_miss30_maf04_noBadInds_noHighCov.recode.vcf
        ##this writes variants_miss30_maf04_noBadInds_noHighCov.recode.mpgl

Calculate coverage

    $ perl /working/lgalland/perl_scripts/coverage_calc_bam.pl variants_miss30_maf04_noBadInds_noHighCov.recode.mpgl PA_ids_188_good_head.txt /working/lgalland/attenuata/bwa/sam_sai/ &
            # this writes out the variants_miss40_noBadInds_noHighCov.recode.mpgl_coverage.csv

Create pntest file

    $ perl /working/lgalland/perl_scripts/gl2genestV1.3.pl variants_miss30_maf04_noBadInds_noHighCov.recode.mpgl mean

## Filter out paralogs

    /working/lgalland/attenuata/bwa/sam_sai/

    $ mkdir HDplot






















