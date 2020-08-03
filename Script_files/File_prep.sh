#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=1,walltime=5:00:00
#PBS -l vmem=16gb
#PBS -M parkanha@iupui.edu
#PBS -m abe
#PBS -N hg19ref
#PBS -j oe

#modules to load
module load bwa/0.7.12
module load samtools/1.9
module load tabix
module load bedtools

#reference, databases and softwares
##Please change the $path according to your sample directory/software path
REF="/$path/Database/GATK/gatk-bundle/hg19/hg19_chr.fa"
PICARD="/$path/picard_2.21.1/picard.jar"
GATK="/$path/gatk-4.1.4.0/gatk"
# https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk
GATK3="/$path/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
#change following based on genome version, download from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
GATK_BUNDLE_DIR="/$path/Database/GATK/gatk-bundle/hg19/"
MILLS=${GATK_BUNDLE_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
PHASE1INDELS=${GATK_BUNDLE_DIR}/1000G_phase1.indels.hg19.sites.vcf
DBSNP=${GATK_BUNDLE_DIR}/dbsnp_138.hg19.vcf

#index reference
#1) for bwa mem
if ! ls ${REF}".bwt" 1> /dev/null 2>&1; then
	${BWA} index -a bwtsw ${REF}
fi
#2) for GATK
REF_BASE="${REF%.*}"
REF_DICT=${REF_BASE}".dict"
if ! ls ${REF_DICT} 1> /dev/null 2>&1; then
	java -jar ${PICARD} CreateSequenceDictionary \
		REFERENCE=${REF} \
		OUTPUT=${REF_DICT}
fi
#3) for GATK
if ! ls ${REF}".fai" 1> /dev/null 2>&1; then
	samtools faidx ${REF}
fi

#####
#Prepare/Create output folders for the analysis
#Make result directory
#Path to result folder location/home directory
cd ~
mkdir output
cd output
mkdir {fastq,qc,bam,result_files,DB_files,sam_index_bam}
cd qc
mkdir {fastqc,multiqc}
cd ../bam
mkdir temp
cd temp
mkdir {bwa,sorted,markdup,bqsr_indelrealign}
cd ../../result_files
mkdir {cnvkit,fpfilter,hsmetrics,manta,strelka_germline,strelka_mutation,vep}
cd ~
######

######
#Prepare bed files for the analysis. Remember to download the bed files to "BED_files" directory
cd BED_files
sort -k1,1V -k2,2n MyelomaPanel1Mutationsv2_final.BED > Mutation.bed
bgzip -c Mutation.bed > Mutation.bed.gz
tabix -f -p bed Mutation.bed.gz

sort -k1,1V -k2,2n MyelomaPanelALL_final.BED > All.bed
bgzip -c All.bed > All.bed.gz
tabix -f -p bed All.bed.gz

#path to picard
cd ~
#BedToIntervalList for HSMetrics
java -jar picard-2.10.0_picard.jar BedToIntervalList I=BED_files/MyelomaPanel2Translocationsv2.BED O=BED_files/Translocation_list.interval_list SD=${GATK_BUNDLE_DIR}/hg19_chr.dict

java -jar picard-2.10.0_picard.jar BedToIntervalList I=BED_files/MyelomaPanel1Mutationsv2_final.BED O=BED_files/Mutation_list.interval_list SD=${GATK_BUNDLE_DIR}/hg19_chr.dict




