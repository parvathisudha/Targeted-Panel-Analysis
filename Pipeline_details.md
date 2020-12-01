# The software used in this pipeline are:
#### Preprocessing: 
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
- [MultiQC](https://multiqc.info/)
- [BWA-MEM](http://bio-bwa.sourceforge.net/)
- [Picard-v2.10.0](https://github.com/broadinstitute/picard/releases/tag/2.10.0)
- [GATK3](https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk) – for indel realignment
- [GATK4](https://github.com/broadinstitute/gatk/releases/tag/4.1.4.0)
 
#### Variant Analysis:
- [Strelka](https://github.com/Illumina/strelka/releases/tag/v2.9.2)
- [MANTA](https://github.com/Illumina/manta/releases/tag/v1.6.0)
- [fpfilter](https://github.com/ckandoth/variant-filter)
- [bam-readcount](https://gist.github.com/ckandoth/87ba44948cb747916f8d#file-build_bam_readcount-txt)
- [VEP-v96](http://grch37.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer)
- [cnvkit](https://cnvkit.readthedocs.io/en/stable/quickstart.html)
 #### For [TarPan Viewer](https://github.com/tcashby/tarpan)
- [SQLite](https://www.sqlite.org/index.html)
- [DB Viewer](https://sqlitebrowser.org/)
- [R](https://www.r-project.org/)
- [R Shiny](https://shiny.rstudio.com/)
- [RStudio](https://www.rstudio.com/) - easiest way to implement using this IDE
- [Python](https://www.python.org/)
 
Reference genome: hg19 - can be downloaded from Broadinstitute’s ftp site:
 ```sh
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/
````
 NOTE: unzip and save to /path_to_Home_directory/Database/GATK/gatk-bundle/hg19 
 Please remember to remove the non-main chromosomes from the genome before preprocessing as MANTA analysis requires bam files with main chromosomes (1-22, X, Y and MT) only. This can be done using “[fasplit](https://bioconda.github.io/recipes/ucsc-fasplit/README.html)”.
 - GATK_BUNDLE
    - Download all the files to the hg19 folder and unzip them.  
```sh
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/
````

### Installation details for Variant analysis software:
##### [Strelka](https://github.com/Illumina/strelka/releases/tag/v2.9.2)
#
```sh
# Download strelka binary
$ wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6x8664.tar.bz2
# Decompress
$ tar xvjf strelka-2.9.2.centos6x8664.tar.bz2
# Run demo to check successful installation
$ bash strelka-2.9.2.centos6x8664/bin/runStrelkaSomaticWorkflowDemo.bash
$ bash strelka-2.9.2.centos6x8664/bin/runStrelkaGermlineWorkflowDemo.bash
```
##### [fpfilter](https://github.com/ckandoth/variant-filter)
- Download the fpfilter.pl script, and view the detailed usage manual:
```sh
$ curl -LO https://github.com/ckandoth/variant-filter/archive/master.zip
$ unzip master.zip
```
##### [bam-readcount](https://gist.github.com/ckandoth/87ba44948cb747916f8d#file-build_bam_readcount-txt)
- Download and unzip the latest commit in the bam-readcount master branch:
- NOTE: In-source builds are not supported, so create a subfolder named build, and build the source in there:
```sh
$ curl -LO https://github.com/genome/bam-readcount/archive/master.tar.gz
$ tar -zxf master.tar.gz
$ mkdir bam-readcount-master/build
$ cd bam-readcount-master/build
$ cmake ../
$ make
```
##### [VEP](http://grch37.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer)
- Follow the instructions of  -    https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html
##### [MANTA](https://github.com/Illumina/manta/releases/tag/v1.6.0)
 - Download and decompress MANTA binary
Refer: https://github.com/Illumina/manta/blob/master/docs/userGuide/installation.md
```sh
$ wget https://github.com/Illumina/manta/releases/tag/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
$ tar -xjf manta-1.6.0.centos6_x86_64.release_src.tar.bz2
$ mkdir build && cd build
$ manta-1.6.0.centos6_x86_64.release_src/configure --jobs=4 
$ make -j4 install
$ manta-1.6.0.centos6_x86_64.release_src/configure --help
````
##### [cnvkit](https://cnvkit.readthedocs.io/en/stable/quickstart.html)
#
````sh
$ conda create -y -p /filepath/conda-env conda package list
$ conda install cnvkit
$ conda create -n cnvkit
$ conda activate cnvkit
$ conda install numpy scipy pandas matplotlib reportlab biopython pyfaidx pysam pyvcf
$ pip install --user -e 
$ cd /path_to_cnvkit/cnvkit/
$ module load r/3.6.0
$ R()        ##Open R
> Install ("Bioconductor::BiocManager")
> library(BiocManager)
> Rscript -e "source('http://callr.org/install#DNAcopy')"
````

## Analysis info:
#### Preprocessing: 
###### Merge Fastq files if they are from different lanes. 
#
##### QC: Check QC for all the sample fastq files.
 - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc): Generates sequence quality metrics [module load fastqc/0.11.5]
 - [MultiQC](https://multiqc.info/): Compile all quality metric logs into single report. 
##### Sequence alignment:
 - [BWA-MEM](http://bio-bwa.sourceforge.net/) [module load bwa/0.7.12]
 - Use BWA-mem for aligning fastq files to human genome assembly hg19. 
##### Sorting the alignment files (BAM), marking duplicates are performed using Picard.
 - [Picard-v2.10.0](https://github.com/broadinstitute/picard/releases/tag/2.10.0):- / path_to_Home_directory /picard-2.10.0_picard.jar
##### Indel realignment and Base quality recalibration using [GATK](https://gatk.broadinstitute.org/hc/en-us).
- NOTE: Indel realignment was part of GATK3 pipeline, but it is not in GATK4 as they incorporated this local realignment of indels into their callers (HaplotypeCaller and Mutect2). If you prefer to use non-GATK caller like Strelka, it is worth it to include it even if it requires this step running within GATK3 version. For BSQR we use GATK4

#### Variant Analysis
##### [Strelka](https://github.com/Illumina/strelka) 
[README.md](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md)
- Strelka calls germline and somatic small variants from mapped sequencing reads. Strelka accepts input read mappings from BAM or CRAM files. It reports all small variant predictions in VCF 4.1 format.
- BED file used for Somatic variant analysis: Mutation.BED
- BED file used for Germline variant analysis: Mutation_translocations.BED
NOTE: BED files need to be sorted and indexed using tabix.
```sh
$ module load tabix
$ module load bedtools
$ sort -k1,1V -k2,2n MyelomaPanel1Mutationsv2_final.BED > Mutation.bed
$ bgzip -c Mutation.bed > Mutation.bed.gz
$ tabix -f -p bed Mutation.bed.gz
$ sort -k1,1V -k2,2n Mutation_Translocations_final.BED > Mutation_Translocations.bed
$ bgzip -c Mutation_Translocations.bed > Mutation_Translocations.bed.gz
$ tabix -f -p bed Mutation_Translocations.bed.gz
## Strelka_Somatic variants – Configuration:
$ strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --tumorBam /$path_to_bam_folder/Sample1_tumor_final.bam --normalBam /$path_to_bam_folder/Sample1_normal_final.bam --referenceFasta /$path_to_hg19_folder/hg19_chr.fa --runDir /$path_to_strelka_mutation_output_folder/Sample1 --exome --callRegions /$path_to_bedfiles_folder/Mutation.bed.gz
# Strelka_Somatic variants – Execution:
$ /$path_to_strelka_mutation_output_folder/Sample1/runWorkflow.py --quiet -m local -j 8
## Filter passed variants and indels from thevcf file:
$ module load vcftools
$ cd /$path_to_strelka_mutation_output_folder/Sample1/results/variants/
$ gunzip *.gz
$ cat somatic.indels.vcf | vcf-annotate -H > somatic.indels_passed.vcf
$ cat somatic.snvs.vcf | vcf-annotate -H > somatic.snvs_passed.vcf
## Strelka_Germline variants – Configuration:
$ strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam /$path_to_bam_folder/Sample1_normal_final.bam --bam /$path_to_bam_folder/Sample1_tumor_final.bam --referenceFasta /$path_to_hg19_folder/hg19_chr.fa --runDir /$path_to_strelka_germline_output_folder/Sample1 --exome --callRegions /$path_to_bedfiles_folder/ bed Mutation_Translocations.bed.gz
## Strelka_Germline variants – Execution:
$ /$path_to_strelka_germline_output_folder/Sample1/runWorkflow.py --quiet -m local -j 8
## filter passed variants
$ cd path_to_strelka_germline_output_folder/Sample1/results/variants/
$ gunzip *.gz
$ cat variants.vcf | vcf-annotate -H variants.vcf > variants_passed.vcf
##Extract AD from vcf files - for "calculating snpdiff.txt"
$ vcftools --vcf variants_passed.vcf --extract-FORMAT-info AD
```
##### [fpfilter](https://github.com/ckandoth/variant-filter)
- FALSE POSITIVE FILTER for variants will improve the precision of variant and mutation calling by removing artifacts associated with short-read alignment. For somatic mutations, generate bam-readcounts with the Tumor BAM. 
> perl variant-filter-master/fpfilter.pl --help
- We used default options, 
```
--var-file          List of variants in VCF, or tab-delimited list of "CHR POS REF VAR"
--readcount-file    The output of bam-readcount for the genomic loci in var-file
--output-file       Output file, tab-delimited list of variants with appended columns for filter status
--min-read-pos      Minimum avg relative distance of variant from start/end of read [$min_read_pos=0.10]
--min-strandedness  Minimum representation of variant allele on each strand [$min_strandedness=0.01]
--min-var-count     Minimum number of variant-supporting reads [$min_var_count=3]
--min-depth         Minimum read depth required across the variant site [$min_depth =8]
--min-var-frac      Minimum variant allele fraction [$min_var_frac=0.05]
--max-mmqs-diff     Maximum difference of mismatch quality sum between var/ref reads (paralog filter) [$max_mmqs_diff=50]
--max-mapqual-diff  Maximum difference of mapping quality between variant and reference reads [$max_mapqual_diff=30]
--max-readlen-diff  Maximum difference of average supporting read length between var/ref reads (paralog filter) [$max_readlen_diff=25]
--min-var-dist-3    Minimum avg distance to effective 3' end of read (real end or Q2) for variant-supporting reads [$min_var_dist_3=0.2]
--max-var-mmqs      Maximum mismatch quality sum of reference-supporting reads [$max_var_mmqs=100]
```
##### [bam-readcount](https://gist.github.com/ckandoth/87ba44948cb747916f8d#file-build_bam_readcount-txt)
- The purpose of this program is to generate metrics at single nucleotide positions.
-  NOTE: bam files needed to be indexed using [samtools](http://www.htslib.org/doc/samtools-index.html).
```sh
$ perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' somatic.passed.snvs.vcf > snvs.var
$ bam-readcount -q1 -b15 -w1 -l snvs.var -f ref/humangenome.fasta sample.bam > snvs.readcount
$ perl variant-filter-master/fpfilter.pl --var-file somatic.passed.snvs.vcf --readcount-file snvs.readcount --output-file snvs.fpfilter
$ awk -v OFS='\t' '{print $1, $2, $8}' snvs.fpfilter > snvs.fpfilter_tab
$ module load tabix
$ bgzip snvs.fpfilter_tab
$ tabix -s 1 -b 2 -e 2 snvs.fpfilter_tab.gz
$ cat somatic.snvs_passed.vcf | vcf-annotate -a snvs.fpfilter_tab.gz -d key=INFO,ID=ANN,Number=1,Type=String,Description='FP filter annotation' -c CHROM,POS,FILTER -H > snvs.fpfilter.vcf
```
##### [VEP](http://grch37.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer)
- Can be performed using online version.
- Select, 
	1) RefSeq transcripts as Transcript database to use.
	2) HGVS identifier
	3) Restrict results to "One selected consequence per variant". 
- Download and save as VCF format.
##### [MANTA](https://github.com/Illumina/manta/releases/tag/v1.6.0)
- MANTA calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs. 
- Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences in VCF 4.1 format.
```sh
## Configuration - example
$ /$path_to_manta_folder/bin/configManta.py --generateEvidenceBam --tumorBam /$path_to_bam_folder/Sample1_tumor_final.bam --normalBam /$path_to_bam_folder/Sample1_normal_final.bam --exome --referenceFasta /$path_to_hg19_folder/hg19_chr.fa --runDir /$path_to_manta_output_folder/Sample1
## Execution
$ /$path_to_manta_output_folder/Sample1/runWorkflow.py 
```
##### [cnvkit](https://cnvkit.readthedocs.io/en/stable/quickstart.html) - Copy number calling pipeline
###### If installing using conda
```sh
$ module unload python
$ module load anaconda/python3.7/2019.03 #[Cahnge the version accordingly]
$ source activate /$path/.conda/envs/cnvkit
$ cd /$path/cnvkit/cnvkit/
$ module load r/3.6.0
$ Rscript -e "source('http://callr.org/install#DNAcopy')"
#Remember to give the annotated BED file here. 
$ /$Path_to_cnvkit_folder/cnvkit/cnvkit.py batch /$path_to_bam_folder/sample1_tumor_final.bam --normal /$path_to_bam_folder/sample1_normal_final.bam --targets /$path_to_bed_files/Mut_Trans_annot.bed --access /$path_to_cnvkit_folder/cnvkit/data/access-5k-mappable.hg19.bed --fasta / $path_to_hg19_folder/hg19_chr.fa --output-reference /$path_to_cnvkit_output_folder /my_reference.cnn --output-dir /$path_to_cnvkit_output_folder/sample1 
```
##### [HSMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard-): using CollectHsMetrics (Picard)
 - This tool requires an aligned SAM or BAM file as well as bait and target interval files in Picard interval_list format. You should use the bait and interval files that correspond to the capture kit that was used to generate the capture libraries for sequencing, which can generally be obtained from the kit manufacturer. If the baits and target intervals are provided in BED format, you can convert them to the Picard interval_list format using Picard's BedToInterval tool. 
 - Refer: 
     - [https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard-]
     - [https://broadinstitute.github.io/picard/picard-metric-definitions.html]
```sh
$ java -jar picard.jar BedToIntervalList \
    I=input.bed \
    O=list.interval_list \
    SD=/path_to_hg19_folder ucsc.hg19.dict

$ java -jar picard.jar CollectHsMetrics \
      I= path_to_bam_folder/input_reads.bam \
      O= path_to_output_folder/output_hs_metrics.txt \
      R=/path_to_hg19_folder/hg19_chr.fa \
      BAIT_INTERVALS=bait.interval_list \
      TARGET_INTERVALS=target.interval_list
```
#### Script files:
- Preprocessing.sh
- Panel_analysis.sh
- snpdiff.R

###### NOTE:  We performed the analysis in the high performance cluster using PBS torque resource manager. 
 - Shell script info:
 ```sh
	#!/bin/bash
	#PBS -M user@iu.edu
	#PBS -l nodes=1:ppn=1,walltime=15:00:00
	#PBS -l vmem=100gb
	#PBS -m abe
	#PBS -N Analysis
	#PBS -j oe
	#PBS -t 1-30
```
#
The jobs were run in parallel for all the samples using PBS -t option. [I run the analysis for 30 samples. Remember to change this according to your sample size.]
The sample names were saved as text file.

### [TarPan](https://github.com/tcashby/tarpan)
 - TarPan Viewer is a tool used to visually inspect targeted panel sequencing data.
```sh
 Open command prompt
$ python -m pip install Pandas
#to upgrade
$ python -m pip install --upgrade pip
```
##### Required files for Database

- BED files:
 	- TargetRegions.BED
 	- Groups.BED
	- Blacklist.BED
 - Variant analysis result files:
    - Mutation File – VEP annotated “somatic.snvs.vcf” and “somatic.indels.vcf” files from Strelka somatic analysis
    - Structural Variant File – “somaticSV.vcf” file from MANTA analysis
    - SNP file  - “snpdiff.txt” generated from Strelka germline analysis using “snpdiff.R”
    - Depth File – “depth.txt” file generated from cnvkit results -  [Generated from normal.targetcoverage.cnn and normal.targetcoverage.cnn]
    - Metrics file - "hsmetrics.txt" files generated using GATK-hsmetrics


###### In command prompt
#
```cs
## change the directory to tarpan folder
$ cd C:\path_to_folder_location\tarpan
$ scripts\create_db.py -db "TarPan.db" -refgen "hg19" -pipeline "Targeted_panel" -targetbed "TargetRegions.BED" -groupbed "Groups.BED" -blacklist "blacklist.BED"
$ scripts\create_db_entry.py -db "TarPan.db" -sampid "Sample1" -normid " Sample1Norm" -mutfile "Samples\ Sample1_snvs.fpfilter_passed_vep.vcf" -muttool "Strelka" -svfile " Samples\ Sample1_somaticSV.vcf" -svtool "Manta" -snp " Samples\ Sample1_snpdiff.txt" -depth " Samples\ Sample1_depth.txt"
```
Note: Repeat the “create_db_entry.py” for all the samples you need for the database. 
##### For [TarPan Viewer](https://github.com/tcashby/tarpan) 
- Open the RStudio
- Clone the repository
- git clone https://github.com/parvathisudha/tarpan.git
- Create a new R Studio project in the cloned directory and open ui.r, server.r, and install.r.
- Install dependencies (this command should work, but you may need to manually intervene)
> source("install.R")

Assuming the dependencies installed correctly, you should now be able to run the application in R Studio.
change the directory to the database directory
>	Run Server.R
