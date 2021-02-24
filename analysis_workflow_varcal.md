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
    dim(coverage) # 543 * 12101 (because first column is ID)

    coverage1 <- coverage[,-1]
    coverage1[1:10, 1:10]
    mean_vect <- vector()
    for (i in 1:543) { mean_vect <- append(mean_vect, mean(as.numeric(coverage1[i,]))) }
    mean(mean_vect)
        ## 8.374944

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

Populations to be included (32 after combinations)

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
    ## YA + YB = YO, Yosemite A + Yosemite B = Yosemite PIAT
    ## CH, China Pines SCI PIMU
    ## CP + RR = CP, Christy Pines + Ridge Road = Christy Pines SCI PIMU
    ## DC, Diablo Canyon PIMU
    ## DM, Del Monte Forest PIMU
    ## FR, Fort Ross PIMU
    ## LO, Lompoc PIMU
    ## NR, Navarro River PIMU
    ## PA, Point Arena PIMU
    ## PB, Pelican Bay SCI PIMU
    ## PP, Partrick Point PIMU
    ## PR, Point Reyes PIMU
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
    grep "PM_CH" pine_543_ids_col.txt > CH_bams.txt
    grep "PM_DC" pine_543_ids_col.txt > DC_bams.txt
    grep "PM_DM" pine_543_ids_col.txt > DM_bams.txt
    grep "PM_FR" pine_543_ids_col.txt > FR_bams.txt
    grep "PM_LO" pine_543_ids_col.txt > LO_bams.txt
    grep "PM_NR" pine_543_ids_col.txt > NR_bams.txt
    grep "PM_PA" pine_543_ids_col.txt > PA_bams.txt
    grep "PM_PB" pine_543_ids_col.txt > PB_bams.txt
    grep "PM_PP" pine_543_ids_col.txt > PP_bams.txt
    grep "PM_PR" pine_543_ids_col.txt > PR_bams.txt
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

    grep "PA_YA" pine_543_ids_col.txt > PA_YA_bams.txt
    grep "PA_YB" pine_543_ids_col.txt > PA_YB_bams.txt

    grep "PM_CP" pine_543_ids_col.txt > PM_CP_bams.txt
    grep "PM_RR" pine_543_ids_col.txt > PM_RR_bams.txt

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
    sed "s/PM_CH/aln_PM_CH/g" CH_bams.txt | sed "s/\$/\.sorted\.bam/g" > CH_bam_names.txt
    sed "s/PM_DC/aln_PM_DC/g" DC_bams.txt | sed "s/\$/\.sorted\.bam/g" > DC_bam_names.txt
    sed "s/PM_DM/aln_PM_DM/g" DM_bams.txt | sed "s/\$/\.sorted\.bam/g" > DM_bam_names.txt
    sed "s/PM_FR/aln_PM_FR/g" FR_bams.txt | sed "s/\$/\.sorted\.bam/g" > FR_bam_names.txt
    sed "s/PM_LO/aln_PM_LO/g" LO_bams.txt | sed "s/\$/\.sorted\.bam/g" > LO_bam_names.txt
    sed "s/PM_NR/aln_PM_NR/g" NR_bams.txt | sed "s/\$/\.sorted\.bam/g" > NR_bam_names.txt
    sed "s/PM_PA/aln_PM_PA/g" PA_bams.txt | sed "s/\$/\.sorted\.bam/g" > PA_bam_names.txt
    sed "s/PM_PB/aln_PM_PB/g" PB_bams.txt | sed "s/\$/\.sorted\.bam/g" > PB_bam_names.txt
    sed "s/PM_PP/aln_PM_PP/g" PP_bams.txt | sed "s/\$/\.sorted\.bam/g" > PP_bam_names.txt
    sed "s/PM_PR/aln_PM_PR/g" PR_bams.txt | sed "s/\$/\.sorted\.bam/g" > PR_bam_names.txt
    sed "s/PM_SP/aln_PM_SP/g" SP_bams.txt | sed "s/\$/\.sorted\.bam/g" > SP_bam_names.txt
    sed "s/PM_SR/aln_PM_SR/g" SR_bams.txt | sed "s/\$/\.sorted\.bam/g" > SR_bam_names.txt
    sed "s/PR_CM/aln_PR_CM/g" CM_bams.txt | sed "s/\$/\.sorted\.bam/g" > CM_bam_names.txt
    sed "s/PR_CN/aln_PR_CN/g" CN_bams.txt | sed "s/\$/\.sorted\.bam/g" > CN_bam_names.txt
    sed "s/PR_CS/aln_PR_CS/g" CS_bams.txt | sed "s/\$/\.sorted\.bam/g" > CS_bam_names.txt
    sed "s/PR_GU/aln_PR_GU/g" GU_bams.txt | sed "s/\$/\.sorted\.bam/g" > GU_bam_names.txt
    sed "s/PR_MP/aln_PR_MP/g" MP_bams.txt | sed "s/\$/\.sorted\.bam/g" > MP_bam_names.txt

combining the proper inds:

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

combining different "populations" that are actually overlapping:

    sed "s/PM_CP/aln_PM_CP/g" PM_CP_bams.txt | sed "s/\$/\.sorted\.bam/g" > PM_CP_bam_names.txt
    sed "s/PM_RR/aln_PM_RR/g" PM_RR_bams.txt | sed "s/\$/\.sorted\.bam/g" > PM_RR_bam_names.txt

    sed "s/PA_YA/aln_PA_YA/g" PA_YA_bams.txt | sed "s/\$/\.sorted\.bam/g" > PA_YA_bam_names.txt
    sed "s/PA_YB/aln_PA_YB/g" PA_YB_bams.txt | sed "s/\$/\.sorted\.bam/g" > PA_YB_bam_names.txt

    cat PM_CP_bam_names.txt PM_RR_bam_names.txt > CP_bam_names.txt
    cat PA_YA_bam_names.txt PA_YB_bam_names.txt > YO_bam_names.txt

    rm PM_RR_bam_names.txt
    rm PM_CP_bam_names.txt
    rm PA_YA_bam_names.txt
    rm PA_YB_bam_names.txt

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
    angsd -bam YO_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out YO &
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
    angsd -bam SP_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out SP &
    angsd -bam SR_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out SR &
    angsd -bam CM_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out CM &
    angsd -bam CS_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out CS &
    angsd -bam GU_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out GU &
    angsd -bam MP_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out MP &
    angsd -bam MR_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out MR &
    angsd -bam SC_bam_names.txt -doSaf 1 -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 -P 2 -out SC &


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
    realSFS CH.saf.idx -P 2 > CH.sfs &
    realSFS DC.saf.idx -P 2 > DC.sfs &
    realSFS DM.saf.idx -P 2 > DM.sfs &
    realSFS FR.saf.idx -P 2 > FR.sfs &
    realSFS LO.saf.idx -P 2 > LO.sfs &
    realSFS NR.saf.idx -P 2 > NR.sfs &
    realSFS PA.saf.idx -P 2 > PA.sfs &
    realSFS PB.saf.idx -P 2 > PB.sfs &
    realSFS PP.saf.idx -P 2 > PP.sfs &
    realSFS PR.saf.idx -P 2 > PR.sfs &
    realSFS SP.saf.idx -P 2 > SP.sfs &
    realSFS SR.saf.idx -P 2 > SR.sfs &
    realSFS CM.saf.idx -P 2 > CM.sfs &
    realSFS CN.saf.idx -P 2 > CN.sfs &
    realSFS CS.saf.idx -P 2 > CS.sfs &
    realSFS GU.saf.idx -P 2 > GU.sfs &
    realSFS MP.saf.idx -P 2 > MP.sfs &  
    realSFS MR.saf.idx -P 2 > MR.sfs &
    realSFS SC.saf.idx -P 2 > SC.sfs &
    realSFS YO.saf.idx -P 2 > YO.sfs &
    realSFS CP.saf.idx -P 2 > CP.sfs &

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
    angsd -bam CH_bam_names.txt -out CH -doThetas 1 -doSaf 1 -pest CH.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam DC_bam_names.txt -out DC -doThetas 1 -doSaf 1 -pest DC.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam DM_bam_names.txt -out DM -doThetas 1 -doSaf 1 -pest DM.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam FR_bam_names.txt -out FR -doThetas 1 -doSaf 1 -pest FR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam LO_bam_names.txt -out LO -doThetas 1 -doSaf 1 -pest LO.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam NR_bam_names.txt -out NR -doThetas 1 -doSaf 1 -pest NR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam PA_bam_names.txt -out PA -doThetas 1 -doSaf 1 -pest PA.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam PB_bam_names.txt -out PB -doThetas 1 -doSaf 1 -pest PB.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam PP_bam_names.txt -out PP -doThetas 1 -doSaf 1 -pest PP.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam PR_bam_names.txt -out PR -doThetas 1 -doSaf 1 -pest PR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam SP_bam_names.txt -out SP -doThetas 1 -doSaf 1 -pest SP.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam SR_bam_names.txt -out SR -doThetas 1 -doSaf 1 -pest SR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam CM_bam_names.txt -out CM -doThetas 1 -doSaf 1 -pest CM.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam CN_bam_names.txt -out CN -doThetas 1 -doSaf 1 -pest CN.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam CS_bam_names.txt -out CS -doThetas 1 -doSaf 1 -pest CS.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam GU_bam_names.txt -out GU -doThetas 1 -doSaf 1 -pest GU.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam MP_bam_names.txt -out MP -doThetas 1 -doSaf 1 -pest MP.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam MR_bam_names.txt -out MR -doThetas 1 -doSaf 1 -pest MR.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam SC_bam_names.txt -out SC -doThetas 1 -doSaf 1 -pest SC.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 & 
    angsd -bam YO_bam_names.txt -out YO -doThetas 1 -doSaf 1 -pest YO.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &
    angsd -bam CP_bam_names.txt -out CP -doThetas 1 -doSaf 1 -pest CP.sfs -anc /working/lgalland/pines_combined/bwa/sam_sai/pine_ref.fasta -GL 1 &


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
    thetaStat do_stat CH.thetas.idx -win 1 -step 1 -outnames CH.thetaWindow
    thetaStat do_stat DC.thetas.idx -win 1 -step 1 -outnames DC.thetaWindow
    thetaStat do_stat DM.thetas.idx -win 1 -step 1 -outnames DM.thetaWindow
    thetaStat do_stat FR.thetas.idx -win 1 -step 1 -outnames FR.thetaWindow
    thetaStat do_stat LO.thetas.idx -win 1 -step 1 -outnames LO.thetaWindow
    thetaStat do_stat NR.thetas.idx -win 1 -step 1 -outnames NR.thetaWindow
    thetaStat do_stat PA.thetas.idx -win 1 -step 1 -outnames PA.thetaWindow
    thetaStat do_stat PB.thetas.idx -win 1 -step 1 -outnames PB.thetaWindow
    thetaStat do_stat PP.thetas.idx -win 1 -step 1 -outnames PP.thetaWindow
    thetaStat do_stat PR.thetas.idx -win 1 -step 1 -outnames PR.thetaWindow
    thetaStat do_stat SP.thetas.idx -win 1 -step 1 -outnames SP.thetaWindow
    thetaStat do_stat SR.thetas.idx -win 1 -step 1 -outnames SR.thetaWindow
    thetaStat do_stat CM.thetas.idx -win 1 -step 1 -outnames CM.thetaWindow
    thetaStat do_stat CN.thetas.idx -win 1 -step 1 -outnames CN.thetaWindow
    thetaStat do_stat CS.thetas.idx -win 1 -step 1 -outnames CS.thetaWindow
    thetaStat do_stat GU.thetas.idx -win 1 -step 1 -outnames GU.thetaWindow
    thetaStat do_stat MP.thetas.idx -win 1 -step 1 -outnames MP.thetaWindow
    thetaStat do_stat MR.thetas.idx -win 1 -step 1 -outnames MR.thetaWindow
    thetaStat do_stat SC.thetas.idx -win 1 -step 1 -outnames SC.thetaWindow
    thetaStat do_stat CP.thetas.idx -win 1 -step 1 -outnames CP.thetaWindow
    thetaStat do_stat YO.thetas.idx -win 1 -step 1 -outnames YO.thetaWindow


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
    thetaStat print CH.thetas.idx > CH.theta_out
    thetaStat print DC.thetas.idx > DC.theta_out
    thetaStat print DM.thetas.idx > DM.theta_out
    thetaStat print FR.thetas.idx > FR.theta_out
    thetaStat print LO.thetas.idx > LO.theta_out
    thetaStat print NR.thetas.idx > NR.theta_out
    thetaStat print PA.thetas.idx > PA.theta_out
    thetaStat print PB.thetas.idx > PB.theta_out
    thetaStat print PP.thetas.idx > PP.theta_out
    thetaStat print PR.thetas.idx > PR.theta_out
    thetaStat print SP.thetas.idx > SP.theta_out
    thetaStat print SR.thetas.idx > SR.theta_out
    thetaStat print CM.thetas.idx > CM.theta_out
    thetaStat print CN.thetas.idx > CN.theta_out
    thetaStat print CS.thetas.idx > CS.theta_out
    thetaStat print GU.thetas.idx > GU.theta_out
    thetaStat print MP.thetas.idx > MP.theta_out
    thetaStat print MR.thetas.idx > MR.theta_out
    thetaStat print SC.thetas.idx > SC.theta_out
    thetaStat print CP.thetas.idx > CP.theta_out
    thetaStat print YO.thetas.idx > YO.theta_out


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

######################################################################

Make the proper pops file for entropy.

    $ cp pine_543_ids_col.txt  pine_ids_543_col.txt
    $ sed -s  "s/[A-Z][A-Z]_*//" pine_ids_543_col.txt > pine_pops_543.txt
        # takes the no-header ID file entries from IDs like PM_DC_0001 to DC_0001

######################################################################	
# RUNNING ENTROPY 
######################################################################

    /working/lgalland/pines_combined/bwa/sam_sai/
    
    $ mkdir entropy
        # we will use this directory later, when we scp the ldak6.txt etc. files back to the node

scp files to laptop to prepare infiles for entropy

    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/pntest_mean_variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.txt /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy
    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/pine_ids_543_col.txt /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy
    $ scp lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/pine_pops_543.txt /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy
  
	
## 1. LDA for starting values

    setwd("/Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy")

    g <- read.table("pntest_mean_variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.txt", header=F)

    dim(g)
    # 12100  543

    names <- read.table("pine_ids_543_col.txt", header=F)
    pops <- read.table("pine_pops_543.txt", header=F)

    nind <- dim(g)[2]
    nloci <- dim(g)[1]

    gmn<-apply(g,1,mean, na.rm=T)
    gmnmat<-matrix(gmn,nrow=nloci,ncol=nind)
    gprime<-g-gmnmat ## remove mean
    gcovarmat<-matrix(NA,nrow=nind,ncol=nind)
    for(i in 1:nind){
    for(j in i:nind){
        if (i==j){
        gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
        }
        else{
        gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
        gcovarmat[j,i]<-gcovarmat[i,j]
        }
    }
    }


    # PCA on the genotype covariance matrix
    pcgcov<-prcomp(x=gcovarmat,center=TRUE,scale=FALSE)
    imp <- summary(pcgcov)
    summary(pcgcov)


    # LDA

    pcgcov->pcg

    library(MASS)

    k2<-kmeans(pcg$x[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    k3<-kmeans(pcg$x[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    k4<-kmeans(pcg$x[,1:5],4,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    k5<-kmeans(pcg$x[,1:5],5,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    k6<-kmeans(pcg$x[,1:5],6,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    k7<-kmeans(pcg$x[,1:5],7,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    k8<-kmeans(pcg$x[,1:5],8,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    k9<-kmeans(pcg$x[,1:5],9,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    k10<-kmeans(pcg$x[,1:5],10,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    k11<-kmeans(pcg$x[,1:5],11,iter.max=10,nstart=10,algorithm="Hartigan-Wong")


    ldak2<-lda(x=pcg$x[,1:5],grouping=k2$cluster,CV=TRUE)
    ldak3<-lda(x=pcg$x[,1:5],grouping=k3$cluster,CV=TRUE)
    ldak4<-lda(x=pcg$x[,1:5],grouping=k4$cluster,CV=TRUE)
    ldak5<-lda(x=pcg$x[,1:5],grouping=k5$cluster,CV=TRUE)
    ldak6<-lda(x=pcg$x[,1:5],grouping=k6$cluster,CV=TRUE)
    ldak7<-lda(x=pcg$x[,1:5],grouping=k7$cluster,CV=TRUE)
    ldak8<-lda(x=pcg$x[,1:5],grouping=k8$cluster,CV=TRUE)
    ldak9<-lda(x=pcg$x[,1:5],grouping=k9$cluster,CV=TRUE)
    ldak10<-lda(x=pcg$x[,1:5],grouping=k10$cluster,CV=TRUE)
    ldak11<-lda(x=pcg$x[,1:5],grouping=k11$cluster,CV=TRUE)


    write.table(round(ldak2$posterior,5),file="ldak2.txt",quote=F,row.names=F,col.names=F)
    write.table(round(ldak3$posterior,5),file="ldak3.txt",quote=F,row.names=F,col.names=F)
    write.table(round(ldak4$posterior,5),file="ldak4.txt",quote=F,row.names=F,col.names=F)
    write.table(round(ldak5$posterior,5),file="ldak5.txt",quote=F,row.names=F,col.names=F)
    write.table(round(ldak6$posterior,5),file="ldak6.txt",quote=F,row.names=F,col.names=F)
    write.table(round(ldak7$posterior,5),file="ldak7.txt",quote=F,row.names=F,col.names=F)
    write.table(round(ldak8$posterior,5),file="ldak8.txt",quote=F,row.names=F,col.names=F)
    write.table(round(ldak9$posterior,5),file="ldak9.txt",quote=F,row.names=F,col.names=F)
    write.table(round(ldak10$posterior,5),file="ldak10.txt",quote=F,row.names=F,col.names=F)
    write.table(round(ldak11$posterior,5),file="ldak11.txt",quote=F,row.names=F,col.names=F)


From command line, obviously:

    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak2.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy
    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak3.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy
    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak4.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy
    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak5.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy
    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak6.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy
    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak7.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy
    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak8.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy
    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak9.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy
    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak10.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy
    $ scp /Users/lanie/lanie/PhD/genomics/pines/combined_allSpecies/entropy/ldak11.txt lgalland@134.197.63.151:/working/lgalland/pines_combined/bwa/sam_sai/entropy


## 2.  Making .mpgl files for entropy

    /working/lgalland/pines_combined/bwa/sam_sai/
    # note that this is NOT in entropy/ ....yet

    $ perl /working/lgalland/perl_scripts/create_entropy_top_2rows.pl pine_ids_543_col.txt 
    $ cat entropy_2rows.txt variants_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.mpgl > pine_redone_entropy.mpgl

Need to add to the top 543 12100 1 ---- First number is number of individuals, second number is number of loci after ALL filtering, third number is just a 1. To do this, just press "enter" at top line, and then enter the 3 values above (entered on the top, newly created line), separated by spaces).

    $ nano pine_redone_entropy.mpgl

## 3. Running entropy

    /working/lgalland/pines_combined/bwa/sam_sai/
    mv pine_redone_entropy.mpgl entropy/

    /working/lgalland/pines_combined/bwa/sam_sai/entropy/

Make subdirectories for each chain and copy all the LDA files (e.g., ldak2.txt, ldak3.txt, etc.) and "redone" mpgl (i.e., pine_redone_entropy.mpgl) file into each of (run0, run1, run2, run3, and run4 directories) - need to run each k group 5 times. Thus, you should start each k2 chain (run0 thru run4), and let them all run concurrently. Then, continue extracting parameters and record DICs for each chain. After all this is done, you will move to k3 entropy. We do this, because you may reach reach outrageous orders of magnitude difference in DIC values, meaning there's no point in continuing all the way through k12

    $ module load entropy/1.2

    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k2.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 2 -q ldak2.txt -m 1 -w 0 &> k2stdout.txt &
    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k3.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 3 -q ldak3.txt -m 1 -w 0 &> k3stdout.txt &
    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k4.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 4 -q ldak4.txt -m 1 -w 0 &> k4stdout.txt &
    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k5.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 5 -q ldak5.txt -m 1 -w 0 &> k5stdout.txt &
    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k6.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 6 -q ldak6.txt -m 1 -w 0 &> k6stdout.txt &
    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k7.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak7.txt -m 1 -w 0 &> k7stdout.txt &
    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k8.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 8 -q ldak8.txt -m 1 -w 0 &> k8stdout.txt &
    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k9.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 9 -q ldak9.txt -m 1 -w 0 &> k9stdout.txt &
    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k10.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 10 -q ldak10.txt -m 1 -w 0 &> k10stdout.txt &
    $ entropy -i pine_redone_entropy.mpgl -o pine_redone_entropy_k11.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 11 -q ldak11.txt -m 1 -w 0 &> k11stdout.txt &

## 4. Extracting parameter estimates

    $ module load entropy/1.2

Extract q estimates

    $ estpost.entropy pine_redone_entropy_k2.hdf5  -p q -s 0 -o q2.txt
    $ estpost.entropy pine_redone_entropy_k3.hdf5  -p q -s 0 -o q3.txt
    $ estpost.entropy pine_redone_entropy_k4.hdf5  -p q -s 0 -o q4.txt
    $ estpost.entropy pine_redone_entropy_k5.hdf5  -p q -s 0 -o q5.txt
    $ estpost.entropy pine_redone_entropy_k6.hdf5  -p q -s 0 -o q6.txt
    $ estpost.entropy pine_redone_entropy_k7.hdf5  -p q -s 0 -o q7.txt
    $ estpost.entropy pine_redone_entropy_k8.hdf5  -p q -s 0 -o q8.txt
    $ estpost.entropy pine_redone_entropy_k9.hdf5  -p q -s 0 -o q9.txt
    $ estpost.entropy pine_redone_entropy_k10.hdf5  -p q -s 0 -o q10.txt











DONE TO HERE!





    










    $ estpost.entropy pine_redone_entropy_k11.hdf5  -p q -s 0 -o q11.txt

Extract gprob estimates from .hdf5 FSmpressed results:

    $ estpost.entropy  pine_redone_entropy_k2.hdf5 -p gprob -s 0 -o gprob2.txt &
    $ estpost.entropy  pine_redone_entropy_k3.hdf5 -p gprob -s 0 -o gprob3.txt &
    $ estpost.entropy  pine_redone_entropy_k4.hdf5 -p gprob -s 0 -o gprob4.txt &
    $ estpost.entropy  pine_redone_entropy_k5.hdf5 -p gprob -s 0 -o gprob5.txt &
    $ estpost.entropy  pine_redone_entropy_k6.hdf5 -p gprob -s 0 -o gprob6.txt &
    $ estpost.entropy  pine_redone_entropy_k7.hdf5 -p gprob -s 0 -o gprob7.txt &
    $ estpost.entropy  pine_redone_entropy_k8.hdf5 -p gprob -s 0 -o gprob8.txt &
    $ estpost.entropy  pine_redone_entropy_k9.hdf5 -p gprob -s 0 -o gprob9.txt &
    $ estpost.entropy  pine_redone_entropy_k10.hdf5 -p gprob -s 0 -o gprob10.txt &




    $ estpost.entropy  pine_redone_entropy_k11.hdf5 -p gprob -s 0 -o gprob11.txt &

Extract DIC estimates from .hdf5 TAmpressed results:

    $ estpost.entropy pine_redone_entropy_k2.hdf5 -s 3 -p deviance
    $ estpost.entropy pine_redone_entropy_k3.hdf5 -s 3 -p deviance
    $ estpost.entropy pine_redone_entropy_k4.hdf5 -s 3 -p deviance
    $ estpost.entropy pine_redone_entropy_k5.hdf5 -s 3 -p deviance
    $ estpost.entropy pine_redone_entropy_k6.hdf5 -s 3 -p deviance
    $ estpost.entropy pine_redone_entropy_k7.hdf5 -s 3 -p deviance
    $ estpost.entropy pine_redone_entropy_k8.hdf5 -s 3 -p deviance
    $ estpost.entropy pine_redone_entropy_k9.hdf5 -s 3 -p deviance
    $ estpost.entropy pine_redone_entropy_k10.hdf5 -s 3 -p deviance



    
    $ estpost.entropy pine_redone_entropy_k11.hdf5 -s 3 -p deviance



















#######################################
## making entropy barplots
#######################################

# First, scp q files (like q2.txt) and an ID file to laptop 
## IMPORTANT: After seeing the first barplot (I guess?), you need to go through and sort the q file and the IDs file (sort in the same way), so that like-colors (ancestral clusters) are together.

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/OM_ids_150_good_noHead.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/entropy/run0/q2.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/entropy/run2/q3.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/entropy/run4/q4.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/entropy/run1/q5.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/entropy/run0/q6.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/entropy/run0/q7.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/entropy/run3/q8.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy

/working/lgalland/rbt/bwa/sam_sai

############ THIS IS FOR AFTER YOU SEE INITIAL BARPLOTS > the following steps rearrange the populations, based on colors we want next together. Should probably have done this manually to save time, but here's the code since I did it. (My ID file is HT_ids_83_col.txt, indicating that it's the ID file with only 83 individuals, which is the number of individuals after all filtering steps. It's in a column format, with no header, like TA_AM_0024, TA_AM_0094, etc. etc.)

cp OM_ids_150_col.txt entropy/

grep "WR_" OM_ids_150_col.txt > WR_inds.txt
grep "BK_" OM_ids_150_col.txt > BK_inds.txt
grep "MD_MD_" OM_ids_150_col.txt > MD_inds.txt
grep "MC_" OM_ids_150_col.txt > MC_inds.txt
grep "TA_" OM_ids_150_col.txt > TA_inds.txt
grep "UT_" OM_ids_150_col.txt > UT_inds.txt
grep "IN_" OM_ids_150_col.txt > IN_inds.txt
grep "TH_" OM_ids_150_col.txt > TH_inds.txt

cat WR_inds.txt BK_inds.txt MD_inds.txt MC_inds.txt TA_inds.txt UT_inds.txt IN_inds.txt TH_inds.txt > names_sorted_for_final_bar.txt
wc -l names_sorted_for_final_bar.txt
#150

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/entropy/names_sorted_for_final_bar.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy

########## Changing the order of the lines in the q file to match the resorted IDs file. Done manually by opening the q2.txt file in text reader and copying and pasting chunks of individuals. Order: AM, ON, TU. SUPER IMPORTANT TO REMEMBER THIS STEP!!!

##############################
# in R to make entropy barplots
##############################

setwd("/Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy")

##############################
## 1) Load data
##############################

names_sort <- read.table("names_sorted_for_final_bar.txt", header=F)
# names_sort <- read.table("OM_ids_150_good_noHead.txt", header=F)

dim(names_sort) # 150 x 1
head(names_sort)

##############################
## 2) Entropy barplots
##############################


######### k2

# run4 had the lowest DIC.
q_dat_sort_2 <- read.csv("q2_run0_sorted_for_bar.txt", header=T)
# q_dat_sort_2 <- read.csv("q2.txt", header=T)

dim(q_dat_sort_2) # 300 x 5
head(q_dat_sort_2)

nind <- 150

cluster21 <- cbind(names_sort, q_dat_sort_2[1:nind,])
cluster22 <- cbind(names_sort, q_dat_sort_2[(1+nind):(nind*2),])

means2 <- data.frame(cbind(cluster21$mean,cluster22$mean))

t.means2 <- t(means2)
##writing table for with better organized data for future plots
ids_means2.txt <- cbind(names_sort, means2)
##write.table(ids_means2.txt, file = "k2_means_ids.txt

quartz(width=8, height=7.5)
par(mar=c(5,5,3,1),mfrow=c(3,1))

colors= c("deepskyblue", "grey30", "coral", "chartreuse1")
barplot(t.means2, col=colors, beside=F, names.arg=names_sort[1:150,], cex.names=.55, las=2, ylab = "Assignment probability",  space=0,lwd=1.5)

######### k3

q_dat_sort_3 <- read.csv("q3_run2_sorted_for_bar.txt", header=T)
dim(q_dat_sort_3) # 450 x 5
head(q_dat_sort_3)

nind <- 150

cluster31 <- cbind(names_sort, q_dat_sort_3[1:nind,])
cluster32 <- cbind(names_sort, q_dat_sort_3[(1+nind):(nind*2),])
cluster33 <- cbind(names_sort, q_dat_sort_3[(1+nind*2):(nind*3),])

means3 <- data.frame(cbind(cluster31$mean,cluster32$mean,cluster33$mean))

t.means3 <- t(means3)
##writing table for with better organized data for future plots
ids_means3.txt <- cbind(names_sort, means3)
##write.table(ids_means3.txt, file = "k3_means_ids.txt", quote=F, row.names=F, col.names=F)

barplot(t.means3, col=colors, beside=F, names.arg=names_sort[1:150,], cex.names=.55, las=2, ylab = "Assignment probability",  space=0,lwd=1.5)

######### k4

q_dat_sort_4 <- read.csv("q4_run4_sorted_for_bar.txt", header=T)
dim(q_dat_sort_4) # 600 x 5
head(q_dat_sort_4)

nind <- 150

cluster41 <- cbind(names_sort, q_dat_sort_4[1:nind,])
cluster42 <- cbind(names_sort, q_dat_sort_4[(1+nind):(nind*2),])
cluster43 <- cbind(names_sort, q_dat_sort_4[(1+nind*2):(nind*3),])
cluster44 <- cbind(names_sort, q_dat_sort_4[(1+nind*3):(nind*4),])

means4 <- data.frame(cbind(cluster43$mean,cluster42$mean, cluster41$mean, cluster44$mean))

t.means4 <- t(means4)
##writing table for with better organized data for future plots
ids_means4.txt <- cbind(names_sort, means4)
#write.table(ids_means4.txt, file = "k4_means_ids.txt", quote=F, row.names=F, col.names=F)

barplot(t.means4, col=colors, beside=F, names.arg=names_sort[1:150,], cex.names=.55, las=2, ylab = "Assignment probability",  space=0,lwd=1.5)

# SAVE QUARTZ WINDOW AS entropy_barplot_3panel.pdf

# copy and paste everything to notes file!
# Nothing but the k2 plot is informative. Thsi is what I expected, but I went forward with all plots through k6 to be thorough in the event we need the plots for supplementary information.


##############################################################
###### PCA on gprobs from entropy. 
##############################################################

/working/lgalland/rbt/bwa/sam_sai/
  
  # 150 lines, with DL DL DL AM AM AM, etc etc

### OM_pops_150.txt is the ID file we need/want for PCA! 

$ sed -s "s/_[A-Z][A-Z]//" OM_ids_150_col.txt > OM_pops3.txt
$ sed -s "s/_[0-9]*//" OM_pops3.txt > OM_pops_150.txt

# Using the gprob2.txt from the run with the lowest DIC (which is gprob2.txt from run0).

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/OM_pops_150.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy/PCA

scp lgalland@ponderosa.biology.unr.edu:/working/lgalland/rbt/bwa/sam_sai/entropy/run0/gprob2.txt /Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy/PCA

##############################################################
## in R, to make entropy gprob PCAs
##############################################################

setwd("/Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy/PCA")

######################
# Read in the data
#####################

# file with ONLY the ID names, in original order matching the gprob file, no header, AM AM AM AM DL DL DL DL etc etc
popIDs <- read.csv("OM_pops_150.txt", header=FALSE)

# gprob file output from entropy, from the run with the lowest k2 DIC
gprobs <- read.csv("gprob2.txt", header=TRUE)
dim(gprobs)
# 150 x 12808 (one more than the actual 12807 loci). the first column is "ind_0" etc etc (which we don't need)
gprobs[1:10,1:10]

gprobs_noname <- gprobs[,-1] # taking off the first column, which is ind_0, etc etc etc
dim(gprobs_noname)
# 150 x 12807. Even though there is technically a 84th row, it isn't counted because we called it a header=TRUE
gprobs_noname[1:10,1:10]

# need to cbind the HT_83_pops.txt and gprob2.txt, so we have the gprob file with IDs for plotting points and colors
gprob <- cbind(popIDs, gprobs_noname)
dim(gprob) # 150 * 12808
gprob[1:10,1:10] # perfect

##############
# run PCA
##############

pcaout <- prcomp(gprobs_noname, center=TRUE, scale=FALSE)
imp <- summary(pcaout)
summary(pcaout)

# Importance of components:
#   PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10    PC11    PC12    PC13    PC14    PC15    PC16    PC17    PC18
# Standard deviation     7.71872 5.58128 5.48321 4.93270 4.82110 4.75564 4.70291 4.62249 4.56722 4.53266 4.45125 4.43147 4.42384 4.40187 4.33745 4.28764 4.26906 4.22419
# Proportion of Variance 0.03059 0.01599 0.01544 0.01249 0.01193 0.01161 0.01136 0.01097 0.01071 0.01055 0.01017 0.01008 0.01005 0.00995 0.00966 0.00944 0.00936 0.00916
# Cumulative Proportion  0.03059 0.04659 0.06202 0.07452 0.08645 0.09807 0.10942 0.12039 0.13110 0.14165 0.15183 0.16191 0.17196 0.18191 0.19157 0.20101 0.21037 0.21953

### IMPORTANT NOTE: If you want to see the PC scores, or investigate the pcaout, do the following:
# to investigate the srtructure:
str(pcaout)
# to view PC scores, which is the "x" in this case:
pcaout$x[,1]

######### PLOTTING PCA ###########
colors <-c(
  "#7c36bb",
  "#d4967d",
  "#c7007f",
  "#d47000",
  "#00702e",
  "#007dd7",
  "#e0c540",
  "#873752")

quartz(width=5, height=10)
par(mar=c(5,5,3,1),mfrow=c(2,1))

#############
# PC 1v2
#############

plot(pcaout$x[,1], pcaout $x[,2], type="n", xlab=paste("PC1 (",(imp$importance[,1][[2]]*100), "% )", sep=""), ylab=paste("PC2 (",(imp$importance[,2][[2]]*100), "% )", sep=""), cex.lab=1.2)

legend("bottomleft", legend=c("Upper Truckee (UT)", "Taylor (TA)", "McKinney (MC)", "Madden (MD)", "Blackwood (BK)", "Ward (WR)", "Incline (IN)", "Third (TH)"), pch=c(16,16,16,16,16,16,16,16), ncol=2, col=colors[1:8], cex=.52)

points(pcaout $x[which(gprob[,1]=="UT"),1], pcaout $x[which(gprob[,1]=="UT"), 2], pch=21, bg=colors[1], cex=1.2)
points(pcaout $x[which(gprob[,1]=="TA"),1], pcaout $x[which(gprob[,1]=="TA"), 2], pch=21, bg=colors[2], cex=1.2)
points(pcaout $x[which(gprob[,1]=="MC"),1], pcaout $x[which(gprob[,1]=="MC"), 2], pch=21, bg=colors[3], cex=1.2)
points(pcaout $x[which(gprob[,1]=="MD"),1], pcaout $x[which(gprob[,1]=="MD"), 2], pch=21, bg=colors[4], cex=1.2)
points(pcaout $x[which(gprob[,1]=="BK"),1], pcaout $x[which(gprob[,1]=="BK"), 2], pch=21, bg=colors[5], cex=1.2)
points(pcaout $x[which(gprob[,1]=="WR"),1], pcaout $x[which(gprob[,1]=="WR"), 2], pch=21, bg=colors[6], cex=1.2)
points(pcaout $x[which(gprob[,1]=="IN"),1], pcaout $x[which(gprob[,1]=="IN"), 2], pch=21, bg=colors[7], cex=1.2)
points(pcaout $x[which(gprob[,1]=="TH"),1], pcaout $x[which(gprob[,1]=="TH"), 2], pch=21, bg=colors[8], cex=1.2)

#############
# PC 1v2, grouped by region
#############
colors2 <-c(
  "#7c36bb",
  "darkgreen",
  "gold")

plot(pcaout$x[,1], pcaout $x[,2], type="n", xlab=paste("PC1 (",(imp$importance[,1][[2]]*100), "% )", sep=""), ylab=paste("PC2 (",(imp$importance[,2][[2]]*100), "% )", sep=""), cex.lab=1.2)

legend("bottomleft", legend=c("South Lake", "West Lake", "North Lake"), pch=c(16,16,16), ncol=1, col=colors2[1:3], cex=.7)

points(pcaout $x[which(gprob[,1]=="UT"),1], pcaout $x[which(gprob[,1]=="UT"), 2], pch=21, bg=colors2[1], cex=1.2)
points(pcaout $x[which(gprob[,1]=="TA"),1], pcaout $x[which(gprob[,1]=="TA"), 2], pch=21, bg=colors2[1], cex=1.2)
points(pcaout $x[which(gprob[,1]=="MC"),1], pcaout $x[which(gprob[,1]=="MC"), 2], pch=21, bg=colors2[2], cex=1.2)
points(pcaout $x[which(gprob[,1]=="MD"),1], pcaout $x[which(gprob[,1]=="MD"), 2], pch=21, bg=colors2[2], cex=1.2)
points(pcaout $x[which(gprob[,1]=="BK"),1], pcaout $x[which(gprob[,1]=="BK"), 2], pch=21, bg=colors2[2], cex=1.2)
points(pcaout $x[which(gprob[,1]=="WR"),1], pcaout $x[which(gprob[,1]=="WR"), 2], pch=21, bg=colors2[2], cex=1.2)
points(pcaout $x[which(gprob[,1]=="IN"),1], pcaout $x[which(gprob[,1]=="IN"), 2], pch=21, bg=colors2[3], cex=1.2)
points(pcaout $x[which(gprob[,1]=="TH"),1], pcaout $x[which(gprob[,1]=="TH"), 2], pch=21, bg=colors2[3], cex=1.2)
 
#save quartz window as PCA_1v2_IndAndGrouped.pdf

#############
# PC 3v4
#############

plot(pcaout$x[,3], pcaout $x[,4], type="n", xlab=paste("PC3 (",(imp$importance[,3][[2]]*100), "% )", sep=""), ylab=paste("PC4 (",(imp$importance[,4][[2]]*100), "% )", sep=""), cex.lab=1.2)
# the 3 and 4 in imp$importance[,3][[2]]*100), "% )", sep=""), ylab=paste("PC4 (",(imp$importance[,4][[2]]*100), refer to the percent of variance that is calculated and printed on the x and y axis. It does not change the placement of the points on the PCA
# to change the points plotted, you have to change each one of the "points" lines below.

legend("bottomleft", legend=c("Blackwood (BK)", "Incline (IN)", "McKinney (MC)", "Madden (MD)", "Taylor (TA)", "Third (TH)", "Upper Truckee (UT)", "Ward (WR)"), pch=c(16,16,16,16,16,16,16,16), ncol=1, col=colors[1:8], cex=.6)

points(pcaout $x[which(gprob[,1]=="UT"),3], pcaout $x[which(gprob[,1]=="UT"), 4], pch=21, bg=colors[1], cex=1.2)
points(pcaout $x[which(gprob[,1]=="TA"),3], pcaout $x[which(gprob[,1]=="TA"), 4], pch=21, bg=colors[2], cex=1.2)
points(pcaout $x[which(gprob[,1]=="MC"),3], pcaout $x[which(gprob[,1]=="MC"), 4], pch=21, bg=colors[3], cex=1.2)
points(pcaout $x[which(gprob[,1]=="MD"),3], pcaout $x[which(gprob[,1]=="MD"), 4], pch=21, bg=colors[4], cex=1.2)
points(pcaout $x[which(gprob[,1]=="BK"),3], pcaout $x[which(gprob[,1]=="BK"), 4], pch=21, bg=colors[5], cex=1.2)
points(pcaout $x[which(gprob[,1]=="WR"),3], pcaout $x[which(gprob[,1]=="WR"), 4], pch=21, bg=colors[6], cex=1.2)
points(pcaout $x[which(gprob[,1]=="IN"),3], pcaout $x[which(gprob[,1]=="IN"), 4], pch=21, bg=colors[7], cex=1.2)
points(pcaout $x[which(gprob[,1]=="TH"),3], pcaout $x[which(gprob[,1]=="TH"), 4], pch=21, bg=colors[8], cex=1.2)


# save the quartz window as rbt_gprob_all.pdf
# copy and paste all R code into notes file!

##########################################
##########################################
# DIVERSITY METRICS FROM ENTROPY GPROBS
##########################################
##########################################

# Calculate expected and observed heterozygosity, and Fis
##########################################

#In R to calculate He, Ho, Fis

setwd("/Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy/diversity_metrics")

##########################################

# manipulate gprob2 file:

gprob_original <- read.delim("gprob2.txt", header=F, sep=",") # standard gprob file from entropy
dim(gprob_original) # 151 * 12808 (because we actually have a header and a sillt first column of ind_0, etc.) 
gprob_original[1:10,1:10]

# rip off the first column and first row (because we do actually have a header):
gprob_noname <- gprob_original[-c(1),-c(1)]
dim(gprob_noname) # 150 * 12807 
gprob_noname[1:10,1:10]

# need gprob file with IDs, so:
pops <- read.delim("OM_pops_150.txt", header=FALSE) # file with no header, only pop abbreviations in a single column, like this: AM AM Am DL DL DL
pops[1:10,] 

gprob0 <- cbind(pops, gprob_noname)
dim(gprob0)
gprob0[1:10,1:10]

ids <- read.delim("OM_ids_150_good_noHead.txt", header = FALSE) # file with single column, no header, inds listed like this: TA_AM_0013 TA_AM_0034 TA_DL_0044 TA_DL_0055
dim(ids)
ids

gprob <- cbind(ids, gprob0)
dim(gprob) # 150 * 12809 (two more than # SNPs)
gprob[1:10,1:10]

pop_order <- read.delim("rbt_pop_list_ordered.txt", header=FALSE) # file with no header, single column, listing populations ONLY ONCE in the order that you want the table output!!!! LIterally, it's this many rows: DL EG UR ON AM KO MO AJ TU (and that's it!)
dim(pop_order)
head(pop_order)


##########################################
## genetic diversity
##########################################

## expected het

npops <- length(pop_order[,1])

het_out <- matrix(0, npops, 3) # the 3 remains a 3 regardles of npops


for (i in 1:npops)
{
  print(pop_order[i,1])
  pop_subset <- subset(gprob_noname, gprob[,2]==as.character(pop_order[i,1]))
  dim_pop_subset <- dim(pop_subset)
  mean_loc <- vector()
  for (j in 1:12807)
  {
    p <- mean(as.numeric(as.character(pop_subset[,j])) / 2 )
    q <- 1 - p
    het <- 1 - (p^2 + q^2)
    mean_loc <- append(mean_loc, het)
  }
  het_out[i,1] <- as.character(pop_order[i,1])
  het_out[i,2] <- dim_pop_subset[1]
  het_out[i,3] <- mean(mean_loc)
  print(i)
}

het_out


## observed het from gprobs (0.9 - 1.1 = HET)

obs_het_out <- matrix(0, npops, 2)

for (i in 1:npops)
{
  pop_sub <- subset(gprob_noname, gprob[,2]==as.character(pop_order[i,1]))
  dim_pop_sub <- dim(pop_sub)
  het_assign <- matrix(0, dim_pop_sub[1], dim_pop_sub[2])
  mean_het_assign <- vector()
  for (j in 1:dim_pop_sub[1])
  {
    for (k in 1:dim_pop_sub[2])
    {
      if (as.numeric(pop_sub[j,k]) > 0.9 && as.numeric(pop_sub[j,k]) < 1.1)	{ het_assign[j,k] <- 1 }
      else											{ het_assign[j,k] <- 0 }
    }
    mean_het_assign <- append(mean_het_assign, mean(het_assign[j,]))
  }
  obs_het_out[i,1] <- as.character(pop_order[i,1])
  obs_het_out[i,2] <- mean(mean_het_assign)
  print(i)
}


## fis = 1 - (Ho / He)

final_div_mat <- matrix(0, 8, 5) # is the second number in parentheses the number of populations? Yes. 10 pops.
for (i in 1:8) # the 1:10 here refers to # populations
{
  fis <- 1 - (as.numeric(obs_het_out[i,2]) / as.numeric(het_out[i,3]))
  final_div_mat[i,1] <- het_out[i,1]
  final_div_mat[i,2] <- het_out[i,2]
  final_div_mat[i,3] <- het_out[i,3]
  final_div_mat[i,4] <- obs_het_out[i,2]
  final_div_mat[i,5] <- fis
}

final_div_mat

write.table(final_div_mat, file="rbt_diversity_out.txt", quote=F, row.names=F, col.names=F)
## this writes out the table with metrics listed in this order: He, Ho, and Fis. Positive Fis indicates heterozygote deficiency.

## Copy and paste everything back into notes file!

##########################################
##########################################

###################### allele frequencies, Fst, Nei's D, NJ tree ####################################


setwd("/Users/lanie/lanie/PhD/genomics/trout/data_analysis/entropy/diversity_metrics")

# read in gprob file from entropy
gprob_original <- read.delim("gprob2.txt", header=T, sep=",") # standard gprob file from entropy, with best DIC
dim(gprob_original) # 150 * 12808 (one more than # of loci, because we actually have a first column of ind_0, etc.) 
gprob_original[1:10,1:10]

# rip off first column
gmat <- gprob_original[,-c(1)]
dim(gmat) # 150 * 12807 ---- yes, because this is number of inds and number of loci
gmat[1:10,1:10]

nind <- 150
nloci <- 12807
gprob <- gmat
gprob[1:10,1:10]

sg.pops_v<-read.table("OM_pops_150.txt", header=F) 
# file with only pops listed in single column, in same order as gprob file. (Example: AM AM AM AM TU TU TU TU TU etc etc to 192 lines)
as.vector(sg.pops_v[,1])->sg.pops

sg.afreq <- matrix(data=NA, ncol=nloci, nrow=length(unique(sg.pops)))

for(i in 1:nloci){
  for(j in 1:length(unique(sg.pops))){
    sg.afreq[j,(i)] <- sum(gprob[(sg.pops == unique(sg.pops)[j]),i])/(2*length(sg.pops[sg.pops == unique(sg.pops)[j]]))
    
  }
}
##afreq matrix below can be used for hudsons fst.

write.table(round(sg.afreq,digits=5), file="rbt_afreq_matrix.txt", quote=F, row.names=F, col.names=F, sep=",")


## This script calculates distance matrix, Nei's D
all.afreq <- read.csv("rbt_afreq_matrix.txt", header=F) # this is the file you'll use to create the nice pairwise table and the histograms of pairwise Fst
## function Nei's D 1983 from Genetic Distances and Reconstruction of Phylogenetic Trees From Microsatellite DNA, Genetics 1996

neisD<-function(P){ ## rows = pops, cols = loci
  N<-dim(P)[1]
  L<-dim(P)[2] 
  D<-matrix(0,nrow=N,ncol=N)
  for(i in 1:(N-1)){ for(j in (i+1):N){
    D[i,j]<-1 - (sum(sqrt((1-P[i,])*(1-P[j,]))) + sum(sqrt(P[i,]*P[j,]))
    )/L
    D[j,i]<-D[i,j]
  }}
  D<-as.dist(D)
  return(D)
}


## calculate nei's distance
neiDcommon<-neisD(all.afreq) # this takes a bit of time. No, it's not frozen!

library("ape")
read.table("rbt_pop_list_ordered.txt", header=F) -> pop_ids # rbt_pop_list_ordered.txt is a file with each population abbreviation listed ONLY ONCE, in a single column. So, you have 9 populations? You should have only 9 rows indicating unique population identifiers
mat.neiD <- as.matrix(neiDcommon)
rownames(mat.neiD) <- pop_ids[,1]
write.table(mat.neiD, file="neisDmatrix_rbt.txt", col.names=F, quote=F)

intermediate <- read.table("neisDmatrix_rbt.txt", header=F)
mat.neiD <- intermediate[,-c(1)]

read.table("rbt_pop_list_ordered.txt", header=F) -> pop_ids
rownames(mat.neiD) <- pop_ids[,1]

library(ape)

tree <- bionj(as.dist(mat.neiD))

pdf("nj_dist_tree.pdf", width=10, height=6)
par(mfrow=c(1,2))

plot(tree, "unrooted", use.edge.length=T, lab4ut="axial", no.margin=T, cex= 0.75)

dev.off()

#########################################################################################
# to calculate Hudson's Fst and the pretty Fst/Nei's D on each diagonal table, and Fst pairwise histograms
#	This code uses a function for hudsons Fst to calculate locus specific Fst   # #estimates for all pairwise combinations of individuals. The input is an allele frequency# #matrix with nrows = the number of individuals and ncols = the number of loci. This code# #will write out a pdf with frequency distribution plots and a file with locus specific # #fst estimates in columns, with a header that is each possible pairwise pop/species     # #combination.	
#########################################################################################

hudsonFst2<-function(p1=NA,p2=NA,n1=NA,n2=NA){
  numerator<-p1 * (1 - p1) + p2 * (1 - p2)
  denominator<-p1 * (1 - p2) + p2 * (1 - p1)
  fst<-1 - numerator/denominator
  out<-cbind(numerator,denominator,fst)
  return(out)
}

species <- c("BK", "IN", "MC", "MD", "TA", "TH", "UT", "WR") ##taxon labels, or in this case, population labels
species <-as.matrix(species)
nloci <- 12807


#MUST SPECIFY THE NUMBER DIMENSIONS OF EXPECTED OUTPUT
pop.fst.mean.med <- matrix(data=NA, nrow=8, ncol=8) #the "10" values here refer to number of populations
af<- read.table("rbt_afreq_matrix.txt", header=F, sep = ",")
afreq<-as.matrix(af)



## Calculate locus-specific Fst, plot histogram, populate matrix of mean/median
results <- data.frame()
names <-matrix()


pdf("loc_fst_hist.pdf", width = 15, height=20)
#quartz(width=2, height=10)
#par(mar=c(5,5,3,1),mfrow=c(8,7))

par(mfrow=c(9,5)) # first number is number of rows
par(mar=c(5,4,4,2))
par(oma=c(2,2,2,2))
for (i in 1:7){ # the 7 listed here is ONE LESS than the number of populations
  for (j in i:8){ # the 8 here is THE SAME AS the number of populations
    fst <- hudsonFst2(p1=as.numeric(afreq[i,]), p2=as.numeric(afreq[j,]))
    if(i != j){
      hist(fst[,3], col="grey", main=paste(toupper(species[i]), toupper(species[j]), sep=","), breaks=seq(0,1,0.025), cex.lab = 1.2, cex = 1.2, cex.main = 1.2, xlab= "Locus-specific Hudson's Fst", ylab="frequency", xlim=c(0,.8))
      abline(v=mean(fst[,3], na.rm=T), lty=2, col="darkred")
      text(0.5, 2500, paste("mean Fst = ", round(mean(fst[,3], na.rm=T), digits=2)))
      box()
      id<-as.character(paste(species[i],species[j],sep="_"))
      
      # mat.fst[]<-fst[,3]
      results <- rbind(results, fst[,3])
      names <- cbind(names, id)
      
    }
    pop.fst.mean.med[i,j] <- mean(fst[,3], na.rm=T)
    pop.fst.mean.med[j,i] <- quantile(fst[,3], probs=c(0.5), na.rm=T)
    
  }
}


names <-t(names)

write.table(results, file = "loc_fst_pairs.csv", row.names = FALSE, 
            append = FALSE, col.names = FALSE, quote=F, sep = ",")
write.table(names, file = "loc_fst_pairs_names.csv", row.names = FALSE, 
            append = FALSE, col.names = FALSE, quote=F, sep = ",")           
dev.off()

colnames(pop.fst.mean.med) <- species
write.csv(cbind(species,pop.fst.mean.med), file="pop_fst_mean_med.csv", col.names = F, quote=F, row.names=F)
# ignore the warning that the previous line yields

# copy and paste all code back to notes file!

####################################################################################


















