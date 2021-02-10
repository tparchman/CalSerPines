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


Calculate the mean number of reads per individual, from your individual fastqs

    grep "^@" -c *.fastq > meanReads_perInd.txt &
    
    R
    meanReads <- read.delim("meanReads_perInd.txt", header=FALSE, sep=":")
    mean(meanReads[,2])
2590771 (this is the mean!)

    quit()

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

`NOTE`: See new_gl_attempt2.R in /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/bwa/PCA for the highly edited PCAs where populations are plotted in helpful order, and individuals are removed.

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

![PCA_pinesCombined](https://github.com/tparchman/CalSerPines/blob/master/pinesCombined_firstLook_miss30_maf05.pdf) 

![PCA_pinesCombined_groupedBySpecies](https://github.com/tparchman/CalSerPines/blob/master/pinesCombined_firstLook_bySpecies.pdf) 

Make coverage file. Give it your  mpgl file, your ids file, and the location of your `bam` files. The output file is *.recode.mpgl_coverage.csv. You can go in another terminal window and verify this file is increasing in size.
	
	perl /working/lgalland/perl_scripts/coverage_calc_bam.pl variants_miss30_maf05_noBadInds.recode.mpgl pine_ids_553_good_head.txt /working/lgalland/pines_combined/bwa/sam_sai/ &
	
Need to filter out over-assembled loci. Start by cutting the scaffold names column from the original file into the new file, for use in R.    

    $ cut -d " " -f 1 variants_miss30_maf05_noBadInds.recode.mpgl > loc_names_15404.txt

Move files for use in R.

    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/variants_miss30_maf05_noBadInds.recode.mpgl_coverage.csv /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies

    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/loc_names_15404.txt /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies

Identify over-assembled loci in R

    R
    
    setwd("/Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies")

    dat <- read.csv("variants_miss30_maf05_noBadInds.recode.mpgl_coverage.csv", header=F)
        dim(dat) #553 * 15405 - yes, because first column is ID
        dat[1:10,1:10]

    loc_names <- read.delim("loc_names_15404.txt", header=F)
        head(loc_names)
        dim(loc_names) #15404 * 1

    dat_noname <- dat[,-1] # rip off the first column, which is ID
        dim(dat_noname) #553 * 15404
        dat_noname[1:10,1:10]

    #RERUN THIS FOLLOWING VECTOR EVERY TIME YOU CHANGE PARAMETERS, BECAUSE THE FOR LOOP BELOW APPENDS AND AMMENDS IT
    #start with 10, 12, 15, 20, 25. If no tail is visible, increase to 25, 30, 35, 40, 45

    avg_15404 <- vector()
    in_out_15404_20 <- vector()
    in_out_15404_25 <- vector()
    in_out_15404_30 <- vector()
    in_out_15404_35 <- vector()
    in_out_15404_40 <- vector()

    for (i in 1:15404)
        {
        avg <- mean(dat_noname[,i])
        avg_15404 <- append(avg_15404, avg)
  
        if (avg <= 20)	{ in_out_15404_20 <- append(in_out_15404_20, 1) }
        else			{ in_out_15404_20 <- append(in_out_15404_20, 0) }
  
        if (avg <= 25)	{ in_out_15404_25 <- append(in_out_15404_25, 1) }
        else			{ in_out_15404_25 <- append(in_out_15404_25, 0) }
  
        if (avg <= 30)	{ in_out_15404_30 <- append(in_out_15404_30, 1) }
        else			{ in_out_15404_30 <- append(in_out_15404_30, 0) }
  
        if (avg <= 35)	{ in_out_15404_35 <- append(in_out_15404_35, 1) }
        else			{ in_out_15404_35 <- append(in_out_15404_35, 0) }
  
        if (avg <= 40)	{ in_out_15404_40 <- append(in_out_15404_40, 1) }
        else			{ in_out_15404_40 <- append(in_out_15404_40, 0) }
        }

    hist(avg_15404)
    hist(avg_15404, xlim=c(0,80), breaks=2000) # can change the scale here
    
    # if you look at this histogram and don't see a tail, change the numbers above (10, 12, 15, 20) to other values in the for loop

    ##choose the value at which the numbers STOP declining.

    sum(in_out_15404_20) 
        ## 12119
    sum(in_out_15404_25) 
        ## 12689
    sum(in_out_15404_30) 
        ## 13086
    sum(in_out_15404_35) 	
        ## 13411
    sum(in_out_15404_40)  ## choosing to kill all locs with mean cov/ind >= 40
        ## 13637


    sub_40 <- dat_noname[,in_out_15404_40==1]
        dim(sub_40)
        # 553 x 13637    
    sub_40_avg <- subset(avg_15404, in_out_15404_40==1)	

    kill_locs <- subset(loc_names, in_out_15404_40==0)
        dim(kill_locs)
        # 1767 x 1 (because 15404 - 1767 = 13637, or the number of loci remaining after we remove loci) 
        head(kill_locs)

    write.table(kill_locs, file="high_cov_loc_list_to_be_removed.txt", quote=F, row.names=F, col.names=F)

    #### copy and paste all of this back to the notes file, to keep record of everything that was done!

Time to filter out these over-assembled loci 

    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/high_cov_loc_list_to_be_removed.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai

    $ sed "s/:/\t/" high_cov_loc_list_to_be_removed.txt > high_cov_loc_list_to_be_removed_tabdelim.txt

    $ module load vcftools/0.1.14 

    $ vcftools --vcf variants_miss30_maf05_noBadInds.recode.vcf --exclude-positions high_cov_loc_list_to_be_removed_tabdelim.txt --recode --recode-INFO-all --out variants_miss30_maf05_noBadInds_noHighCov
	    ## After filtering, kept 13637 out of a possible 15404 Sites

Generate NEW depth of coverage files, with no bad inds, no over-assembled loci
    $ vcftools --vcf variants_miss30_maf05_noBadInds_noHighCov.recode.vcf --depth -c > variants_miss30_maf05_noBadInds_noHighCov.txt

Create mpgl file NEW .mpgl file (no bad inds, no over-assembled loci)

    $ perl /working/lgalland/perl_scripts/vcf2mpglV1.3TLP.pl variants_miss30_maf05_noBadInds_noHighCov.recode.vcf
        ##this writes variants_miss30_maf05_noBadInds_noHighCov.recode.mpgl

Calculate coverage

    $ perl /working/lgalland/perl_scripts/coverage_calc_bam.pl variants_miss30_maf05_noBadInds_noHighCov.recode.mpgl pine_ids_553_good_head.txt /working/lgalland/pines_combined/bwa/sam_sai/ &
            # this writes out the variants_miss40_noBadInds_noHighCov.recode.mpgl_coverage.csv
 
Create pntest file

    $ perl /working/lgalland/perl_scripts/gl2genestV1.3.pl variants_miss30_maf05_noBadInds_noHighCov.recode.mpgl mean

## Filter out paralogs

    /working/lgalland/pines_combined/bwa/sam_sai/

    $ mkdir HDplot

Most updated vcf:

    $ cp variants_miss30_maf05_noBadInds_noHighCov.recode.vcf HDplot/

Copy relevant scripts

    $ cp /working/lgalland/python_scripts/HDplot_python.py /working/lgalland/pines_combined/bwa/sam_sai/HDplot
    $ cp /working/lgalland/python_scripts/vcf_to_depth.py /working/lgalland/pines_combined/bwa/sam_sai/HDplot

    HDplot/

Alter the python script to call the correct [vcf_file =] and [depths_file =]
    
    $ nano HDplot_python.py

Run HDplot

    $ source activate py27
    $ python HDplot_python.py
            ### this will begin rapidly printing numbers to the screen. Finishes quickly. Spits out HDplot_python_exampleInput.depthsBias, HDplot_results.pdf, H_vs_AlleleRatio_results.pdf, variants_miss30_maf04_noBadInds_noHighCov.recode.depths, and vcf_to_depth.pyc. We only really care about H_vs_AlleleRatio_results.pdf and HDplot_results.pdf.

    $ source deactivate

scp these three files to laptop

    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/HDplot/HDplot_python_exampleInput.depthsBias /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/HDplot/

    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/HDplot/HDplot_results.pdf /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/HDplot/

    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/HDplot/H_vs_AlleleRatio_results.pdf /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/HDplot/
        # Look at the pdfs (especially HDplot_results.pdf), and choose cutoffs. Here, we are going to call 0.55 our cutoff for proportion of heterozygotes and we will retain loci between 35 and -30 (D, or read ratio deviation)

In the terminal, working in R, to subset singleton loci based on cutoffs outlined above - from Katie

    /working/lgalland/pines_combined/bwa/sam_sai/HDplot/
    
    R

    HD <- read.table("HDplot_python_exampleInput.depthsBias", header = TRUE) ## 13636
    dim(HD) ## 13636
    #change cutoffs as necessary:
    made_z_cutoff <- intersect(which(HD$z <= 35), which(HD$z >= -30))
    made_z_and_h_cutoff <- intersect(which(HD$hetPerc < 0.55), made_z_cutoff) 

    singleton_snps <- data.frame(HD$contig[made_z_and_h_cutoff], HD$pos[made_z_and_h_cutoff])

    write.table(singleton_snps, file = "singleton_snps.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    quit()

Exiting R, to actually remove singletons! (from Katie)

    /working/lgalland/pines_combined/bwa/sam_sai/HDplot/

    $ mv singleton_snps.txt /working/lgalland/pines_combined/bwa/sam_sai/
    $ cd /working/lgalland/pines_combined/bwa/sam_sai/

    $ less singleton_snps.txt
    $ wc -l singleton_snps.txt
        # 12100 snps remaining. Started with 13636, so we removed 1536 snps!!! Yikes, 11.3% of snps to remove 

Filter paralogs with vcftools

    $ module load vcftools/0.1.14

    $ vcftools --vcf variants_miss30_maf05_noBadInds_noHighCov.recode.vcf --positions singleton_snps.txt --recode --out variants_miss30_maf05_noBadInds_noHighCov_noParalogs
	    ## After filtering, kept 12100 out of a possible 13637 Sites

        ## name of new file: variants_miss30_maf05_noBadInds_noHighCov_noParalogs.recode.vcf

NEED TO REMOVE WEIRD SAMPLES THAT DON'T MAKE SENSE (mislabeled, slightly contaminated, perhaps?). Then, I'll calculate diversity, coverage, etc. Looked at various PCAs and determined the following samples needed to be removed. 

    $ touch weirdOnes.indv
    $ nano weirdOnes.indv

Grab the IDs from out.imiss and copy and paste into weirdOnes.indv

        INDV
        PA_BS_0139,553953
        PA_SL_0119,801470
        PA_SL_0120,757491
        PM_PB_0281,2585951
        PM_PP_0153,502374
        PM_PR_0172,410688
        PM_CP_0340,901310
        PM_DC_0216,793819
        PR_CN_0749,488616
        PR_MP_0827,1378882

Remove the weird ones.

    $ vcftools --vcf variants_miss30_maf05_noBadInds_noHighCov_noParalogs.recode.vcf --remove weirdOnes.indv --recode --recode-INFO-all --out variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird
			After filtering, kept 543 out of 553 Individuals
            After filtering, kept 12100 out of a possible 12100 Sites

Need to make new IDs and pops files that do not include these weird individuals. Easiest to do this manually. Remove the following individuals from EACH of the ID and pop files below.
    
    PA_BS_0139,553953
    PA_SL_0119,801470
    PA_SL_0120,757491
    PM_PB_0281,2585951
    PM_PP_0153,502374
    PM_PR_0172,410688
    PM_CP_0340,901310
    PM_DC_0216,793819
    PR_CN_0749,488616
    PR_MP_0827,1378882

    $ cp pine_553_ids_col.txt pine_543_ids_col.txt
    $ nano pine_543_ids_col.txt
    
    $ cp pine_ids_553_good_head.txt pine_ids_543_good_head.txt
    $ nano pine_ids_543_good_head.txt

    $ cp pine_ids_553_good_noHead_original.txt pine_ids_543_good_noHead_original.txt
    $ nano pine_ids_543_good_noHead_original.txt

    $ cp pine_ids_553_good_noHead.txt pine_ids_543_good_noHead.txt
    $ nano pine_ids_543_good_noHead.txt

    $ cp pine_ids_553.txt pine_ids_543.txt
    $ nano pine_ids_543.txt

Then, `sed` out the parts you don't want for your pops file (removing the "first term" in your ID file, in this case). (The following worked here, which was in format XX(species)_XX(population)_XXXX(ind).) 

	$ sed -s "s/[A-Z][A-Z]_//" pine_543_ids_col.txt > pine_pops3.txt
	$ sed -s "s/_[0-9]*//" pine_pops3.txt > pine_543_pops.txt


Generate NEW depth of coverage, mpgl, and pntest files (no bad inds, no over-assembled loci, no paralogs, no weird ones)

    $ vcftools --vcf variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.vcf --depth -c > variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.txt

    $ perl /working/lgalland/perl_scripts/vcf2mpglV1.3TLP.pl variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.vcf
        ##this writes .mpgl file

    $ module load bwa/0.7.5a
    $ module load samtools/1.3

    $ perl /working/lgalland/perl_scripts/coverage_calc_bam.pl variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.mpgl pine_ids_543_good_head.txt /working/lgalland/pines_combined/bwa/sam_sai/ &
        # creates *.recode.mpgl_coverage.csv
 
    $ perl /working/lgalland/perl_scripts/gl2genestV1.3.pl variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.mpgl mean
        # creates pntest file

Calculate mean coverage. Read final filtered coverage file into R - mean coverage per individual

    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.mpgl_coverage.csv /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/

In R

    R

    setwd("/Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/")

    coverage <- read.csv("variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.mpgl_coverage.csv", header=F)
    coverage[1:10,1:10]
    dim(coverage) # 543 * 8563 (because first column is ID)

    coverage1 <- coverage[,-1]
    coverage1[1:10, 1:10]
    mean_vect <- vector()
    for (i in 1:543) { mean_vect <- append(mean_vect, mean(as.numeric(coverage1[i,]))) }
    mean(mean_vect)
        ## 12.28467

################################################################################################

### diversity with angsd (trevor figured all this out initially) 


Relevant citations:

	angsd:			
        Korneliussen, T. S., Albrechtsen, A., & Nielsen, R. (2014). ANGSD: analysis of next generation sequencing data. BMC Bioinformatics, 15, 356.
	D:				
        Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.
	pi:				
        Tajima, F. (1983). Evolutionary relationship of DNA sequences in finite populations. Genetics, 105(2), 437-460.
	theta w:			
        Watterson, G. A. (1975). On the number of segregating sites in genetical models without recombination. Theoretical Population Biology, 7, 256-276.
	this pipeline:	
        Korneliussen, T. S., Moltke, I., Albrechtsen, A., & Nielsen, R. (2013). Calculation of Tajima’s D and other neutrality test statistics from low depth next-generation sequencing data. BMC Bioinformatics, 14, 289.


Create conda for angsd (Korneliussen et al. 2014 BMC Bioinformatics) (if you don't already have it)

    $ conda create -n angsd python=3.4 anaconda
    $ conda install -c bioconda angsd
    $ conda install --override-channels -c conda-forge bzip2

Get rolling:

    $ source activate angsd
    $ source deactivate angsd

Populations to be included
	## AH, Auburn High PIMU			
	## AL, Aubur Low PIMU
    ## BS, Big Sur PIAT
    ## CG, Cuesta Grade PIAT
    ## LA, Los Angeles PIAT
    ## LS, Lake Shasta PIAT
    ## OC, Orange County PIAT
    ## PF, Panther Flat PIAT
    ## SB, San Bernardino PIAT
    ## SH, Santa Cruz High PIAT
    ## ST, Santa Cruz Top PIAT
    ## YA, Yosemite A PIAT
    ## YB, Yosemite B PIAT
    ## CH, China Pines SCI PIMU
    ## CP, Christy Pines SCI PIMU
    ## DC, Diablo Canyon PIMU
    ## DM, Del Monte Forest PIMU
    ## FR, Fort Ross PIMU
    ## LO, Lompoc PIMU
    ## NR, Navarro River PIMU
    ## PA, Point Arena PIMU
    ## PB, Pelican Bay SCI PIMU
    ## PP, Partrick Point PIMU
    ## PR, Point Reyes PIMU
    ## RR, Ridge Road SCI PIMU
    ## SP, Salt Point PIMU
    ## SR, Sea Ranch PIMU
    ## CM, Cambria PIRA
    ## CN, Cedros North PIRA
    ## CS, Cedros South PIRA
    ## GU, Isla Guadalupe PIRA
    ## MP, Monterey Peninsula PIRA
    ## MR, Monterey PIRA
    ## SC, Santa Cruz PIRA

Make bam lists per pop

    grep "PA_AH" pine_543_ids_col.txt > AH_bams.txt
    grep "PA_AL" pine_543_ids_col.txt > AL_bams.txt
    grep "PA_BS" pine_543_ids_col.txt > BS_bams.txt
    grep "PA_CG" pine_543_ids_col.txt > CG_bams.txt
    grep "PA_LA" pine_543_ids_col.txt > LA_bams.txt
    grep "PA_LS" pine_543_ids_col.txt > LS_bams.txt
    grep "PA_OC" pine_543_ids_col.txt > OC_bams.txt
    grep "PA_PF" pine_543_ids_col.txt > PF_bams.txt
    grep "PA_SB" pine_543_ids_col.txt > SB_bams.txt
    grep "PA_SH" pine_543_ids_col.txt > SH_bams.txt
    grep "PA_ST" pine_543_ids_col.txt > ST_bams.txt
    grep "PA_YA" pine_543_ids_col.txt > YA_bams.txt
    grep "PA_YB" pine_543_ids_col.txt > YB_bams.txt
    grep "PM_CH" pine_543_ids_col.txt > CH_bams.txt
    grep "PM_CP" pine_543_ids_col.txt > CP_bams.txt
    grep "PM_DC" pine_543_ids_col.txt > DC_bams.txt
    grep "PM_DM" pine_543_ids_col.txt > DM_bams.txt
    grep "PM_FR" pine_543_ids_col.txt > FR_bams.txt
    grep "PM_LO" pine_543_ids_col.txt > LO_bams.txt
    grep "PM_NR" pine_543_ids_col.txt > NR_bams.txt
    grep "PM_PA" pine_543_ids_col.txt > PA_bams.txt
    grep "PM_PB" pine_543_ids_col.txt > PB_bams.txt
    grep "PM_PP" pine_543_ids_col.txt > PP_bams.txt
    grep "PM_PR" pine_543_ids_col.txt > PR_bams.txt
    grep "PM_RR" pine_543_ids_col.txt > RR_bams.txt
    grep "PM_SP" pine_543_ids_col.txt > SP_bams.txt
    grep "PM_SR" pine_543_ids_col.txt > SR_bams.txt
    grep "PR_CM" pine_543_ids_col.txt > CM_bams.txt
    grep "PR_CN" pine_543_ids_col.txt > CN_bams.txt
    grep "PR_CS" pine_543_ids_col.txt > CS_bams.txt
    grep "PR_GU" pine_543_ids_col.txt > GU_bams.txt
    grep "PR_MP" pine_543_ids_col.txt > MP_bams.txt

    grep "PR_MR" pine_543_ids_col.txt > PR_MR_bams.txt
    grep "PX_MR" pine_543_ids_col.txt > PX_MR_bams.txt

    grep "PR_SC" pine_543_ids_col.txt > PR_SC_bams.txt
    grep "PX_SC" pine_543_ids_col.txt > PX_SC_bams.txt


    sed "s/PA_AH/aln_PA_AH/g" AH_bams.txt | sed "s/\$/\.sorted\.bam/g" > AH_bam_names.txt
    sed "s/PA_AL/aln_PA_AL/g" AL_bams.txt | sed "s/\$/\.sorted\.bam/g" > AL_bam_names.txt
    sed "s/PA_BS/aln_PA_BS/g" BS_bams.txt | sed "s/\$/\.sorted\.bam/g" > BS_bam_names.txt
    sed "s/PA_CG/aln_PA_CG/g" CG_bams.txt | sed "s/\$/\.sorted\.bam/g" > CG_bam_names.txt
    sed "s/PA_LA/aln_PA_LA/g" LA_bams.txt | sed "s/\$/\.sorted\.bam/g" > LA_bam_names.txt
    sed "s/PA_LS/aln_PA_LS/g" LS_bams.txt | sed "s/\$/\.sorted\.bam/g" > LS_bam_names.txt
    sed "s/PA_OC/aln_PA_OC/g" OC_bams.txt | sed "s/\$/\.sorted\.bam/g" > OC_bam_names.txt
    sed "s/PA_PF/aln_PA_PF/g" PF_bams.txt | sed "s/\$/\.sorted\.bam/g" > PF_bam_names.txt
    sed "s/PA_SB/aln_PA_SB/g" SB_bams.txt | sed "s/\$/\.sorted\.bam/g" > SB_bam_names.txt
    sed "s/PA_SH/aln_PA_SH/g" SH_bams.txt | sed "s/\$/\.sorted\.bam/g" > SH_bam_names.txt
    sed "s/PA_ST/aln_PA_ST/g" ST_bams.txt | sed "s/\$/\.sorted\.bam/g" > ST_bam_names.txt
    sed "s/PA_YA/aln_PA_YA/g" YA_bams.txt | sed "s/\$/\.sorted\.bam/g" > YA_bam_names.txt
    sed "s/PA_YB/aln_PA_YB/g" YB_bams.txt | sed "s/\$/\.sorted\.bam/g" > YB_bam_names.txt
    sed "s/PM_CH/aln_PM_CH/g" CH_bams.txt | sed "s/\$/\.sorted\.bam/g" > CH_bam_names.txt
    sed "s/PM_CP/aln_PM_CP/g" CP_bams.txt | sed "s/\$/\.sorted\.bam/g" > CP_bam_names.txt
    sed "s/PM_DC/aln_PM_DC/g" DC_bams.txt | sed "s/\$/\.sorted\.bam/g" > DC_bam_names.txt
    sed "s/PM_DM/aln_PM_DM/g" DM_bams.txt | sed "s/\$/\.sorted\.bam/g" > DM_bam_names.txt
    sed "s/PM_FR/aln_PM_FR/g" FR_bams.txt | sed "s/\$/\.sorted\.bam/g" > FR_bam_names.txt
    sed "s/PM_LO/aln_PM_LO/g" LO_bams.txt | sed "s/\$/\.sorted\.bam/g" > LO_bam_names.txt
    sed "s/PM_NR/aln_PM_NR/g" NR_bams.txt | sed "s/\$/\.sorted\.bam/g" > NR_bam_names.txt
    sed "s/PM_PA/aln_PM_PA/g" PA_bams.txt | sed "s/\$/\.sorted\.bam/g" > PA_bam_names.txt
    sed "s/PM_PB/aln_PM_PB/g" PB_bams.txt | sed "s/\$/\.sorted\.bam/g" > PB_bam_names.txt
    sed "s/PM_PP/aln_PM_PP/g" PP_bams.txt | sed "s/\$/\.sorted\.bam/g" > PP_bam_names.txt
    sed "s/PM_PR/aln_PM_PR/g" PR_bams.txt | sed "s/\$/\.sorted\.bam/g" > PR_bam_names.txt
    sed "s/PM_RR/aln_PM_RR/g" RR_bams.txt | sed "s/\$/\.sorted\.bam/g" > RR_bam_names.txt
    sed "s/PM_SP/aln_PM_SP/g" SP_bams.txt | sed "s/\$/\.sorted\.bam/g" > SP_bam_names.txt
    sed "s/PM_SR/aln_PM_SR/g" SR_bams.txt | sed "s/\$/\.sorted\.bam/g" > SR_bam_names.txt
    sed "s/PR_CM/aln_PR_CM/g" CM_bams.txt | sed "s/\$/\.sorted\.bam/g" > CM_bam_names.txt
    sed "s/PR_CN/aln_PR_CN/g" CN_bams.txt | sed "s/\$/\.sorted\.bam/g" > CN_bam_names.txt
    sed "s/PR_CS/aln_PR_CS/g" CS_bams.txt | sed "s/\$/\.sorted\.bam/g" > CS_bam_names.txt
    sed "s/PR_GU/aln_PR_GU/g" GU_bams.txt | sed "s/\$/\.sorted\.bam/g" > GU_bam_names.txt
    sed "s/PR_MP/aln_PR_MP/g" MP_bams.txt | sed "s/\$/\.sorted\.bam/g" > MP_bam_names.txt

    sed "s/PR_MR/aln_PR_MR/g" PR_MR_bams.txt | sed "s/\$/\.sorted\.bam/g" > PR_MR_bam_names.txt
    sed "s/PX_MR/aln_PX_MR/g" PX_MR_bams.txt | sed "s/\$/\.sorted\.bam/g" > PX_MR_bam_names.txt
    cat PR_MR_bam_names.txt PX_MR_bam_names.txt > MR_bam_names.txt

    sed "s/PX_SC/aln_PX_SC/g" PX_SC_bams.txt | sed "s/\$/\.sorted\.bam/g" > PX_SC_bam_names.txt
    sed "s/PR_SC/aln_PR_SC/g" PR_SC_bams.txt | sed "s/\$/\.sorted\.bam/g" > PR_SC_bam_names.txt
    cat PR_SC_bam_names.txt PX_SC_bam_names.txt > SC_bam_names.txt

    rm PR_MR_bam_names.txt
    rm PX_MR_bam_names.txt
    rm PX_SC_bam_names.txt
    rm PR_SC_bam_names.txt


## First step (doSaf)
	## -bam <INPUT> = input of bam names for population
	## -doSaf <INPUT> = option 1: calculate the site allele frequency likelihood based on individual genotype likelihoods assuming HWE
	## -anc <INPUT> = ancestral fasta file (i.e. the genome)
	## -GL <INPUT> = genotype likelihood model (1 = SAMtools; 2 = GATK; 3 = SOAPsnp; 4 = SYK)
	## -P <INPUT> = number of cores used
	## -out <INPUT> = outfile prefix
	
#testing -P cores here. NOTE: Stick with -P 2!!!!! If you go higher, everything goes WAY slower.

    angsd -bam AH_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out AH &
    angsd -bam AL_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out AL &
    angsd -bam BS_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out BS &
    angsd -bam CG_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out CG &
    angsd -bam LA_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out LA &
    angsd -bam LS_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out LS &
    angsd -bam OC_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out OC &
    angsd -bam PF_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out PF &
    angsd -bam SB_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out SB &
    angsd -bam SH_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out SH &
    angsd -bam ST_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out ST &
    angsd -bam YA_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out YA &
    angsd -bam YB_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out YB &
    angsd -bam CH_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out CH &
    angsd -bam CP_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out CP &
    angsd -bam DC_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out DC &
    angsd -bam DM_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out DM &
    angsd -bam FR_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out FR &
    angsd -bam LO_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out LO &   
    angsd -bam NR_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out NR &
    angsd -bam PA_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out PA &
    angsd -bam PB_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out PB &
    angsd -bam PP_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out PP &
    angsd -bam PR_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out PR &
    angsd -bam RR_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out RR &
    angsd -bam SP_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out SP &
    angsd -bam SR_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out SR &
    angsd -bam CM_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out CM &
    angsd -bam CS_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out CS &
    angsd -bam GU_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out GU &
    angsd -bam MP_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out MP &
    angsd -bam MR_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out MR &
    angsd -bam SC_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out SC &









# DONE TO HERE










## Second step (realSFS)

	## estimates the (multi) SFS based on a .saf.idx file generated from step 1 (doSaf)
	## -P <INPUT> = number of cores used

    realSFS AH.saf.idx -P 2 > AH.sfs &
    realSFS AL.saf.idx -P 2 > AL.sfs &
    realSFS BS.saf.idx -P 2 > BS.sfs &
    realSFS CG.saf.idx -P 2 > CG.sfs &
    realSFS LA.saf.idx -P 2 > LA.sfs &
    realSFS LS.saf.idx -P 2 > LS.sfs &
    realSFS OC.saf.idx -P 2 > OC.sfs &
    realSFS PF.saf.idx -P 2 > PF.sfs &
    realSFS SB.saf.idx -P 2 > SB.sfs &
    realSFS SH.saf.idx -P 2 > SH.sfs &
    realSFS ST.saf.idx -P 2 > ST.sfs &
    realSFS YA.saf.idx -P 2 > YA.sfs &
    realSFS YB.saf.idx -P 2 > YB.sfs &
    realSFS CH.saf.idx -P 2 > CH.sfs &
    realSFS CP.saf.idx -P 2 > CP.sfs &
    realSFS DC.saf.idx -P 2 > DC.sfs &
    realSFS DM.saf.idx -P 2 > DM.sfs &
    realSFS FR.saf.idx -P 2 > FR.sfs &
    realSFS LO.saf.idx -P 2 > LO.sfs &
    realSFS NR.saf.idx -P 2 > NR.sfs &
    realSFS PA.saf.idx -P 2 > PA.sfs &
    realSFS PB.saf.idx -P 2 > PB.sfs &
    realSFS PP.saf.idx -P 2 > PP.sfs &
    realSFS PR.saf.idx -P 2 > PR.sfs &
    realSFS RR.saf.idx -P 2 > RR.sfs &
    realSFS SP.saf.idx -P 2 > SP.sfs &
    realSFS SR.saf.idx -P 2 > SR.sfs &
    realSFS CM.saf.idx -P 2 > CM.sfs &
    realSFS CN.saf.idx -P 2 > CN.sfs &
    realSFS CS.saf.idx -P 2 > CS.sfs &
    realSFS GU.saf.idx -P 2 > GU.sfs &
    realSFS MP.saf.idx -P 2 > MP.sfs &  
    realSFS MR.saf.idx -P 2 > MR.sfs &
    realSFS SC.saf.idx -P 2 > SC.sfs &

## Third step (doThetas)

	## calculate the thetas (population scale mutation rates) for each site
	## -bam <INPUT> = input of bam names for population
	## -out <INPUT> = outfile prefix
	## -doThetas 1 = calculate the thetas (use 1, not sure what other option values exist or do)
	## -doSaf <INPUT> = option 1: calculate the site allele frequency likelihood based on individual genotype likelihoods assuming HWE
	## -pest = if -pest is supplied in addition to -doSaf then the output will then be posterior probability of the sample allele frequency for each site
	## -anc <INPUT> = ancestral fasta file (i.e. the genome)
	## -GL <INPUT> = genotype likelihood model (1 = SAMtools; 2 = GATK; 3 = SOAPsnp; 4 = SYK)

    cd /working/lgalland/pines_combined/bwa/sam_sai/
    source activate angsd

    angsd -bam AH_bam_names.txt -out AH -doThetas 1 -doSaf 1 -pest AH.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam AL_bam_names.txt -out AL -doThetas 1 -doSaf 1 -pest AL.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam BS_bam_names.txt -out BS -doThetas 1 -doSaf 1 -pest BS.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam CG_bam_names.txt -out CG -doThetas 1 -doSaf 1 -pest CG.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam LA_bam_names.txt -out LA -doThetas 1 -doSaf 1 -pest LA.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam LS_bam_names.txt -out LS -doThetas 1 -doSaf 1 -pest LS.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam OC_bam_names.txt -out OC -doThetas 1 -doSaf 1 -pest OC.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam PF_bam_names.txt -out PF -doThetas 1 -doSaf 1 -pest PF.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam SB_bam_names.txt -out SB -doThetas 1 -doSaf 1 -pest SB.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam SH_bam_names.txt -out SH -doThetas 1 -doSaf 1 -pest SH.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam ST_bam_names.txt -out ST -doThetas 1 -doSaf 1 -pest ST.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam YA_bam_names.txt -out YA -doThetas 1 -doSaf 1 -pest YA.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam YB_bam_names.txt -out YB -doThetas 1 -doSaf 1 -pest YB.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam CH_bam_names.txt -out CH -doThetas 1 -doSaf 1 -pest CH.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam CP_bam_names.txt -out CP -doThetas 1 -doSaf 1 -pest CP.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam DC_bam_names.txt -out DC -doThetas 1 -doSaf 1 -pest DC.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam DM_bam_names.txt -out DM -doThetas 1 -doSaf 1 -pest DM.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam FR_bam_names.txt -out FR -doThetas 1 -doSaf 1 -pest FR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam LO_bam_names.txt -out LO -doThetas 1 -doSaf 1 -pest LO.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam NR_bam_names.txt -out NR -doThetas 1 -doSaf 1 -pest NR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam PA_bam_names.txt -out PA -doThetas 1 -doSaf 1 -pest PA.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam PB_bam_names.txt -out PB -doThetas 1 -doSaf 1 -pest PB.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam PP_bam_names.txt -out PP -doThetas 1 -doSaf 1 -pest PP.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam PR_bam_names.txt -out PR -doThetas 1 -doSaf 1 -pest PR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam RR_bam_names.txt -out RR -doThetas 1 -doSaf 1 -pest RR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam SP_bam_names.txt -out SP -doThetas 1 -doSaf 1 -pest SP.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam SR_bam_names.txt -out SR -doThetas 1 -doSaf 1 -pest SR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam CM_bam_names.txt -out CM -doThetas 1 -doSaf 1 -pest CM.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam CN_bam_names.txt -out CN -doThetas 1 -doSaf 1 -pest CN.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam CS_bam_names.txt -out CS -doThetas 1 -doSaf 1 -pest CS.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam GU_bam_names.txt -out GU -doThetas 1 -doSaf 1 -pest GU.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam MP_bam_names.txt -out MP -doThetas 1 -doSaf 1 -pest MP.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam MR_bam_names.txt -out MR -doThetas 1 -doSaf 1 -pest MR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam SC_bam_names.txt -out SC -doThetas 1 -doSaf 1 -pest SC.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 


## Fourth step (thetaStat do_stat)
	## estimate Tajimas D and other statistics
	## -win = window size
	## -step = step size (step must be equal or greater than win if you want to avoid overlap)
	## D stored in *thetaWindow.pestPG
	## IMPORTANT: BASED ON YOUR REFERENCE, make a decision about window and step size. If using a reference genome, set to -win 50000 -step 50000. If it's a de novo assembly, use -win 1 -step 1. This is because our contigs are....

    thetaStat do_stat AH.thetas.idx -win 1 -step 1 -outnames AH.thetaWindow
    thetaStat do_stat AL.thetas.idx -win 1 -step 1 -outnames AL.thetaWindow
    thetaStat do_stat BS.thetas.idx -win 1 -step 1 -outnames BS.thetaWindow
    thetaStat do_stat CG.thetas.idx -win 1 -step 1 -outnames CG.thetaWindow
    thetaStat do_stat LA.thetas.idx -win 1 -step 1 -outnames LA.thetaWindow
    thetaStat do_stat LS.thetas.idx -win 1 -step 1 -outnames LS.thetaWindow
    thetaStat do_stat OC.thetas.idx -win 1 -step 1 -outnames OC.thetaWindow
    thetaStat do_stat PF.thetas.idx -win 1 -step 1 -outnames PF.thetaWindow
    thetaStat do_stat SB.thetas.idx -win 1 -step 1 -outnames SB.thetaWindow
    thetaStat do_stat SH.thetas.idx -win 1 -step 1 -outnames SH.thetaWindow
    thetaStat do_stat ST.thetas.idx -win 1 -step 1 -outnames ST.thetaWindow
    thetaStat do_stat YA.thetas.idx -win 1 -step 1 -outnames YA.thetaWindow
    thetaStat do_stat YB.thetas.idx -win 1 -step 1 -outnames YB.thetaWindow
    thetaStat do_stat CH.thetas.idx -win 1 -step 1 -outnames CH.thetaWindow
    thetaStat do_stat CP.thetas.idx -win 1 -step 1 -outnames CP.thetaWindow
    thetaStat do_stat DC.thetas.idx -win 1 -step 1 -outnames DC.thetaWindow
    thetaStat do_stat DM.thetas.idx -win 1 -step 1 -outnames DM.thetaWindow
    thetaStat do_stat FR.thetas.idx -win 1 -step 1 -outnames FR.thetaWindow
    thetaStat do_stat LO.thetas.idx -win 1 -step 1 -outnames LO.thetaWindow
    thetaStat do_stat NR.thetas.idx -win 1 -step 1 -outnames NR.thetaWindow
    thetaStat do_stat PA.thetas.idx -win 1 -step 1 -outnames PA.thetaWindow
    thetaStat do_stat PB.thetas.idx -win 1 -step 1 -outnames PB.thetaWindow
    thetaStat do_stat PP.thetas.idx -win 1 -step 1 -outnames PP.thetaWindow
    thetaStat do_stat PR.thetas.idx -win 1 -step 1 -outnames PR.thetaWindow
    thetaStat do_stat RR.thetas.idx -win 1 -step 1 -outnames RR.thetaWindow
    thetaStat do_stat SP.thetas.idx -win 1 -step 1 -outnames SP.thetaWindow
    thetaStat do_stat SR.thetas.idx -win 1 -step 1 -outnames SR.thetaWindow
    thetaStat do_stat CM.thetas.idx -win 1 -step 1 -outnames CM.thetaWindow
    thetaStat do_stat CN.thetas.idx -win 1 -step 1 -outnames CN.thetaWindow
    thetaStat do_stat CS.thetas.idx -win 1 -step 1 -outnames CS.thetaWindow
    thetaStat do_stat GU.thetas.idx -win 1 -step 1 -outnames GU.thetaWindow
    thetaStat do_stat MP.thetas.idx -win 1 -step 1 -outnames MP.thetaWindow
    thetaStat do_stat MR.thetas.idx -win 1 -step 1 -outnames MR.thetaWindow
    thetaStat do_stat SC.thetas.idx -win 1 -step 1 -outnames SC.thetaWindow



## Fourth step data extraction (Tajima's D)
	
    ## Based on script that trevor made
	## Done in R on ponderosa

    R

    library(data.table)

    CI <- function(data){
        low <- mean(data) - 1.96*(sd(data)/sqrt(length(data)))
        high <- mean(data) + 1.96*(sd(data)/sqrt(length(data)))
        return(c(low,high))
    }

# reads in all files in 
    
    pestPG_files <- list.files(pattern='*Window.pestPG')
  
# makes empty data from to fill

    TajD_df <- as.data.frame(matrix(nrow=length(pestPG_files),ncol=4))
    names(TajD_df) <- c('Pop','TajD','TajD_low','TajD_high')
    
    for (i in 1:length(pestPG_files)){
        infile <- fread(pestPG_files[i])
        Pop <- unlist(strsplit(pestPG_files[i],'.',fixed=TRUE))[1]
        TajD_df$Pop[i] <- Pop
        
        #TajD
        TajD <- infile$Tajima
        TajD_df$TajD[i] <- mean(TajD)  
        TajD_ci <- CI(TajD)
        TajD_df$TajD_low[i] <- TajD_ci[1]
        TajD_df$TajD_high[i] <- TajD_ci[2]
        
    }
    rm('infile')
    write.csv(TajD_df,'full_angsd_TajD_out.csv',row.names=FALSE)

    quit()


## Fifth step (thetaStat print)
	## chromosome
	## position
	## Theta Watterson (theta w, Watterson 1975 = Number of segregating sites)
	## Theta Nucleotide diversity (theta pi, Tajima 1983 = average pairwise differences among individuals)
	## Theta Singleton category
	## Theta H
	## Theta L

    thetaStat print AH.thetas.idx > AH.theta_out
    thetaStat print AL.thetas.idx > AL.theta_out
    thetaStat print BS.thetas.idx > BS.theta_out
    thetaStat print CG.thetas.idx > CG.theta_out
    thetaStat print LA.thetas.idx > LA.theta_out
    thetaStat print LS.thetas.idx > LS.theta_out
    thetaStat print OC.thetas.idx > OC.theta_out
    thetaStat print PF.thetas.idx > PF.theta_out
    thetaStat print SB.thetas.idx > SB.theta_out
    thetaStat print SH.thetas.idx > SH.theta_out
    thetaStat print ST.thetas.idx > ST.theta_out
    thetaStat print YA.thetas.idx > YA.theta_out
    thetaStat print YB.thetas.idx > YB.theta_out
    thetaStat print CH.thetas.idx > CH.theta_out
    thetaStat print CP.thetas.idx > CP.theta_out
    thetaStat print DC.thetas.idx > DC.theta_out
    thetaStat print DM.thetas.idx > DM.theta_out
    thetaStat print FR.thetas.idx > FR.theta_out
    thetaStat print LO.thetas.idx > LO.theta_out
    thetaStat print NR.thetas.idx > NR.theta_out
    thetaStat print PA.thetas.idx > PA.theta_out
    thetaStat print PB.thetas.idx > PB.theta_out
    thetaStat print PP.thetas.idx > PP.theta_out
    thetaStat print PR.thetas.idx > PR.theta_out
    thetaStat print RR.thetas.idx > RR.theta_out
    thetaStat print SP.thetas.idx > SP.theta_out
    thetaStat print SR.thetas.idx > SR.theta_out
    thetaStat print CM.thetas.idx > CM.theta_out
    thetaStat print CN.thetas.idx > CN.theta_out
    thetaStat print CS.thetas.idx > CS.theta_out
    thetaStat print GU.thetas.idx > GU.theta_out
    thetaStat print MP.thetas.idx > MP.theta_out
    thetaStat print MR.thetas.idx > MR.theta_out
    thetaStat print SC.thetas.idx > SC.theta_out



## Fifth step data extraction (pi and theta W)
	
    ## Based on script that trevor made
	## Done in R on ponderosa

    R

    library(data.table)

    CI <- function(data){
    low <- mean(data) - 1.96*(sd(data)/sqrt(length(data)))
    high <- mean(data) + 1.96*(sd(data)/sqrt(length(data)))
    return(c(low,high))
    }

    subsetTheta_files <- list.files(pattern='*theta_out')

    subsetPisum_df <- as.data.frame(matrix(nrow=length(subsetTheta_files),ncol=7))
    names(subsetPisum_df) <- c('Pop','Pi','Pi_low','Pi_high', 'Watt', 'Watt_low', 'Watt_high')


    subsetPisamp_df <- as.data.frame(matrix(ncol=2))
    names(subsetPisamp_df) <- c('Pop','Pi')
    for (i in 1:length(subsetTheta_files)){
    infile <- fread(subsetTheta_files[i])
    Pop <- unlist(strsplit(subsetTheta_files[i],'.',fixed=TRUE))[1]
    pi <- exp(infile$Pairwise) #in log form so must exp
    watt <- exp(infile$Watterson)
    
    subsetPisum_df$Pop[i] <- Pop
    subsetPisum_df$Pi[i] <- mean(pi) #in log form 
    subsetPisum_df$Watt[i] <- mean(watt) #in log form 
    
    ci_pi <- CI(pi)
    subsetPisum_df$Pi_low[i] <- ci_pi[1]
    subsetPisum_df$Pi_high[i] <- ci_pi[2]
    
    ci_watt <- CI(watt)
    subsetPisum_df$Watt_low[i] <- ci_watt[1]
    subsetPisum_df$Watt_high[i] <- ci_watt[2]
    }
    
    rm('infile')
    write.csv(subsetPisum_df,'full_angsd_piWatt_out.csv',row.names=FALSE)

    quit()

###### scp files to laptop to view more easily:

    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/full_angsd_piWatt_out.csv /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/pi_angsd/

    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/full_angsd_TajD_out.csv /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/pi_angsd/




















