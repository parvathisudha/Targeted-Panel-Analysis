#!/bin/bash
#PBS -M parkanha@iu.edu
#PBS -l nodes=1:ppn=1,walltime=15:00:00
#PBS -l vmem=50gb
#PBS -m abe
#PBS -N align
#PBS -j oe
#PBS -t 1-30

module load bwa/0.7.12
module load samtools/1.9
module load fastqc/0.11.5

#folders
#Please change the $path according to your sample directory/software path.
RAW_FASTQ_DIR="/$path/fastq/"
RESULTS_DIR="/$path/output"
#If original fastqs are not merged then please uncomment the next line and change the BWA-mem command accordingly.
#FASTQ_DIR=${RESULTS_DIR}"/fastq/"
FASTQC_DIR=${RESULTS_DIR}"/qc/fastqc/"
MULTIQC_DIR=${RESULTS_DIR}"/qc/multiqc/"
BAM_DIR=${RESULTS_DIR}"/bam/"
BAM_BWA_DIR=${RESULTS_DIR}"/bam/temp/bwa/"
BAM_SORTED_DIR=${RESULTS_DIR}"/bam/temp/sorted/"
BAM_DUPMARK_DIR=${RESULTS_DIR}"/bam/temp/markdup/"
BAM_GATK_DIR=${RESULTS_DIR}"/bam/temp/bqsr_indelrealign/"

#master table
MASTER_SAMPLE="/$path/samples_list.txt"
#reference, databases and softwares
REF="/Database/GATK/gatk-bundle/hg19/hg19_chr.fa"
PICARD="/$path/picard-2.10.0_picard.jar"
GATK="/$path/gatk-4.1.4.0/gatk"
# https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk
GATK3="/$path/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
#change following based on genome version, download from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
GATK_BUNDLE_DIR="/$path/GATK/gatk-bundle/hg19/"
MILLS=${GATK_BUNDLE_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
PHASE1INDELS=${GATK_BUNDLE_DIR}/1000G_phase1.indels.hg19.sites.vcf
DBSNP=${GATK_BUNDLE_DIR}/dbsnp_138.hg19.vcf

sample=$(sed -n ${PBS_ARRAYID}p sample_list.txt)
echo $sample
SAMPLE_ID=$(echo "${sample}" | cut -f 2)
BATCH_ID=$(echo "${sample}" | cut -f 4)
PLATFORM=$(echo "${sample}" | cut -f 5)
CAPTURE=$(echo "${sample}" | cut -f 6)
SEQUENCER=$(echo "${sample}" | cut -f 7)
	
	#If original fastqs are not merged then
	#FASTQ_R1=${FASTQ_DIR}${SAMPLE_ID}"_R1.fastq.gz"
	#FASTQ_R2=${FASTQ_DIR}${SAMPLE_ID}"_R2.fastq.gz"
	
	FASTQ_R1=${RAW_FASTQ_DIR}${SAMPLE_ID}"_R1.fastq.gz"
	FASTQ_R2=${RAW_FASTQ_DIR}${SAMPLE_ID}"_R2.fastq.gz"
	BAM_BWA=${BAM_BWA_DIR}${SAMPLE_ID}"_bwa.bam"
	BAM_SORTED=${BAM_SORTED_DIR}${SAMPLE_ID}"_sorted.bam"
	BAM_DUPMARK=${BAM_DUPMARK_DIR}${SAMPLE_ID}"_dupmarked.bam"
	BAM_DUPMARK_METRICS=${BAM_DUPMARK_DIR}${SAMPLE_ID}"_dupmarked_metrics.txt"
	BAM_REALIGN_TARGETS=${BAM_GATK_DIR}${SAMPLE_ID}"_realignment_targets.list"
	BAM_REALIGN=${BAM_GATK_DIR}${SAMPLE_ID}"_realigned.bam"
	BAM_BQSR_RECALTABLE=${BAM_GATK_DIR}${SAMPLE_ID}"_recal_data.table"
	BAM_FINAL=${BAM_DIR}${SAMPLE_ID}"_final.bam"
	BAM_RMDUP=${BAM_DIR}${SAMPLE_ID}"_final_remdup.bam"
	
	#If original fastqs are not merged then - please uncomment the following.
	#for R in "R1" "R2"; do
	# prepare, compress and run fastqc
	 #   if ! ls ${FASTQ_DIR}${SAMPLE_ID}"_"${R}".fastq.gz" 1> /dev/null 2>&1; then
	 #	    echo ${SAMPLE_ID}" - "${R}" - merging fastq files:"
	 #	    touch ${FASTQ_DIR}${SAMPLE_ID}"_"${R}".fastq.gz"
	 #	    for rawfile in ${RAW_FASTQ_DIR}${SAMPLE_ID}*"_"${R}"_"*; do
	 #	        echo ${rawfile}
	 #	        cat ${rawfile} >> ${FASTQ_DIR}${SAMPLE_ID}"_"${R}".fastq.gz"
	 #	    done
	 #	    echo ${SAMPLE_ID}" - "${R}" - fastqc"
	 #	    ${FASTQC} ${FASTQ_DIR}${SAMPLE_ID}"_"${R}".fastq.gz" -o ${FASTQC_DIR}
	 #	fi
	# done
		    
	
	#run fastqc
	for R in "R1" "R2"; do
	echo ${SAMPLE_ID}" - "${R}" - fastqc"
		    fastqc ${RAW_FASTQ_DIR}${SAMPLE_ID}"_"${R}".fastq.gz" -o ${FASTQC_DIR}
		done
	
	#align to human genome reference (BWA MEM)
	#mkdir -p ${BAM_BWA_DIR}
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_BWA} 1> /dev/null 2>&1; then
		echo ${SAMPLE_ID}" - mapping to human genome reference"
		RG='@RG\tID:'${BATCH_ID}'\tSM:'${SAMPLE_ID}'\tPL:${PLATFORM}\tLB:'${CAPTURE}'\tPU:'${SEQUENCER}
		bwa mem \
	       -t 5 \
	        -M \
	        -R ${RG} \
	        ${REF} \
	        ${FASTQ_R1} \
	        ${FASTQ_R2} | samtools view -bS - > ${BAM_BWA}
	fi
	fi

	#sort alignments (PICARD)
	#mkdir -p ${BAM_SORTED_DIR}
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_SORTED} 1> /dev/null 2>&1; then
		echo ${SAMPLE_ID}" - sorting bam"
		java -Xmx20g -jar ${PICARD} SortSam \
			INPUT=${BAM_BWA} \
			OUTPUT=${BAM_SORTED} \
			SORT_ORDER=coordinate
	fi
	fi

	#mark duplicates (PICARD)
	#mkdir -p ${BAM_DUPMARK_DIR}
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_DUPMARK} 1> /dev/null 2>&1; then
		echo ${SAMPLE_ID}" - marking duplicates in bam"
		java -Xmx20g -jar ${PICARD} MarkDuplicates \
			INPUT=${BAM_SORTED} \
			OUTPUT=${BAM_DUPMARK} \
			METRICS_FILE=${BAM_DUPMARK_METRICS}
		#index bam (PICARD)
		java -Xmx20g -jar ${PICARD} BuildBamIndex \
			INPUT=${BAM_DUPMARK}
	fi
	fi

	#indel realignment
	#mkdir -p ${BAM_GATK_DIR}
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_REALIGN} 1> /dev/null 2>&1; then
		#find targets
		echo ${SAMPLE_ID}" - finding targets for indel realignment"
		java -Xmx25g -jar \
		    ${GATK3} \
		        -T RealignerTargetCreator \
		        -R ${REF} \
		        -known ${MILLS} \
		        -known ${PHASE1INDELS} \
		        -I ${BAM_DUPMARK} \
		        -o ${BAM_REALIGN_TARGETS}
		#apply realignment
		echo ${SAMPLE_ID}" - indel realignment application"
		java -Xmx25g -jar \
		    ${GATK3} \
		        -T IndelRealigner \
		        -R ${REF} \
		        -known ${MILLS} \
		        -known ${PHASE1INDELS} \
		        -LOD 0.4 \
		        --maxReadsForRealignment 10000000 \
		        --maxConsensuses 300 \
		        --maxReadsForConsensuses 1200 \
		        -targetIntervals ${BAM_REALIGN_TARGETS} \
		        -I ${BAM_DUPMARK} \
		        -o ${BAM_REALIGN}
	fi
	fi
	
	#base quality scores recalibration (BQSR, GATK)
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
		#build BQSR model
		echo ${SAMPLE_ID}" - bqsr model"
	    ${GATK} BaseRecalibrator \
	        -R ${REF} \
	        --known-sites ${MILLS} \
	        --known-sites ${PHASE1INDELS} \
	        --known-sites ${DBSNP} \
	        -I ${BAM_REALIGN} \
	        -O ${BAM_BQSR_RECALTABLE}
		#apply BQSR
		echo ${SAMPLE_ID}" - bqsr application"
	    ${GATK} ApplyBQSR \
	        -R ${REF} \
	        -bqsr ${BAM_BQSR_RECALTABLE} \
	        -I ${BAM_REALIGN} \
	        -O ${BAM_FINAL}
	fi
