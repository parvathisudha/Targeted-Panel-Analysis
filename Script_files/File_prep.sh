#!/bin/bash
#PBS -k o
#PBS -l nodes=2:ppn=6,walltime=5:00:00
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
REF="//Database/GATK/gatk-boundle/hg19/ucsc.hg19.fasta"
PICARD="/N/u/parkanha/Carbonate/picard_2.21.1/picard.jar"
GATK="/N/u/parkanha/Carbonate/gatk-4.1.4.0/gatk"
# https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk
GATK3="/N/u/parkanha/Carbonate/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
#change following based on genome version, download from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
GATK_BUNDLE_DIR="/N/u/parkanha/Carbonate/Database/GATK/gatk-boundle/hg19/"
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
cd ~
mkdir output
cd output
mkdir qc
cd qc
mkdir fastqc
mkdir multiqc
cd ..
mkdir bam
cd bam
mkdir temp
cd temp
mkdir bwa
mkdir sorted
mkdir markdup
mkdir bqsr_indelrealign
cd ~

######

######
#Preppare bed files for the analysis
cd BED_files
sort -k1,1V -k2,2n MyelomaPanel1Mutationsv2_final.BED > Mutation.bed
bgzip -c Mutation.bed > Mutation.bed.gz
tabix -f -p bed Mutation.bed.gz

sort -k1,1V -k2,2n MyelomaPanelALL_final.BED > All.bed
bgzip -c All.bed.gz.bed > All.bed.gz
tabix -f -p bed All.bed.gz
cd ~

#BedToIntervalList
java -jar picard-2.10.0_picard.jar BedToIntervalList I=BED_files/MyelomaPanel1Mutationsv2_final.BED O=bedfiles/Translocation_list.interval_list SD=${REF}/ucsc.hg19.dict

java -jar picard-2.10.0_picard.jar BedToIntervalList I=BED_files/Mutation_v2.bed O=bedfiles/Mutation_list.interval_list SD=${REF}/ucsc.hg19.dict



