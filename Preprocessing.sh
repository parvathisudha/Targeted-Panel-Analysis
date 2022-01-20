#!/bin/bash
#SBATCH --mail-user=youremail
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=15:00:00
#SBATCH --mem=30G
#SBATCH -J Test
#SBATCH -o Test_%j.txt
#SBATCH -e Test_%j.err
#SBATCH --array=1-10

#modules to load
module load bwa/0.7.12
module load samtools/1.9
module load fastqc/0.11.5
module load java

#load the config file
source ./Panel_config.sh
echo $DIR
echo $REF
echo $PICARD
echo $GATK_BUNDLE_DIR
echo $FASTQ_DIR

MILLS=$GATK_BUNDLE_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf
PHASE1INDELS=$GATK_BUNDLE_DIR/1000G_phase1.snps.high_confidence.hg38.vcf
DBSNP=$GATK_BUNDLE_DIR/Homo_sapiens_assembly38.dbsnp138.vcf

#Path to output directories
FASTQC_DIR=$DIR/qc/fastqc
BAM_DIR=$DIR/bam
BAM_BWA_DIR=$DIR/bam/temp/bwa
BAM_SORTED_DIR=$DIR/bam/temp/sorted
BAM_DUPMARK_DIR=$DIR/bam/temp/markdup
BAM_GATK_DIR=$DIR/bam/temp/bqsr_indelrealign
HSMETRICS_DIR=$DIR/result_files/hsmetrics

sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt)
echo "Sample_info:-"$sample
FILE_ID=$(echo "${sample}" | cut -f 2)
SAMPLE_ID=$(echo "${sample}" | cut -f 3)
BATCH_ID=$(echo "${sample}" | cut -f 5)
PLATFORM=$(echo "${sample}" | cut -f 6)
CAPTURE=$(echo "${sample}" | cut -f 7)
SEQUENCER=$(echo "${sample}" | cut -f 8)

FASTQ_R1=$FASTQ_DIR/$FILE_ID/$SAMPLE_ID"_R1_001.fastq.gz"
FASTQ_R2=$FASTQ_DIR/$FILE_ID/$SAMPLE_ID"_R2_001.fastq.gz"
BAM_BWA=$BAM_BWA_DIR/$SAMPLE_ID"_bwa.bam"
BAM_SORTED=$BAM_SORTED_DIR/$SAMPLE_ID"_sorted.bam"
BAM_DUPMARK=$BAM_DUPMARK_DIR/$SAMPLE_ID"_dupmarked.bam"
BAM_DUPMARK_METRICS=$BAM_DUPMARK_DIR/$SAMPLE_ID"_dupmarked_metrics.txt"
BAM_REALIGN_TARGETS=$BAM_GATK_DIR/$SAMPLE_ID"_realignment_targets.list"
BAM_REALIGN=$BAM_GATK_DIR/$SAMPLE_ID"_realigned.bam"
BAM_BQSR_RECALTABLE=$BAM_GATK_DIR/$SAMPLE_ID"_recal_data.table"
BAM_FINAL=$BAM_DIR/$SAMPLE_ID"_final.bam"
BAM_RMDUP=$BAM_DIR/$SAMPLE_ID"_final_remdup.bam"

echo ${FASTQ_R1}
echo ${FASTQ_R2}

#Check the QC of the FASTQ files by running fastqc
for R in "R1" "R2"; do
	echo ${SAMPLE_ID}" - "${R}" - fastqc"
	fastqc $FASTQ_DIR/${FILE_ID}/${SAMPLE_ID}"_"${R}"_001.fastq.gz" -o ${FASTQC_DIR}
	done

#Align the paired FASTQ to human genome reference using BWA MEM
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_BWA} 1> /dev/null 2>&1; then
	echo ${SAMPLE_ID}" - mapping to human genome reference"
	RG='@RG\tID:'${BATCH_ID}'\tSM:'${SAMPLE_ID}'\tPL:'${PLATFORM}'\tLB:'${CAPTURE}'\tPU:'${SEQUENCER}
	echo ${RG}
		bwa mem \
		-t 5 \
		-M \
		-R ${RG} \
		${REF} \
		${FASTQ_R1} \
		${FASTQ_R2} | samtools view -bS - > ${BAM_BWA}
	fi
	fi

#Sort aligned BAMs using PICARD
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_SORTED} 1> /dev/null 2>&1; then
	echo ${SAMPLE_ID}" - sorting bam"
	java -Xmx20g -jar ${PICARD} SortSam \
		INPUT=${BAM_BWA} \
		OUTPUT=${BAM_SORTED} \
		SORT_ORDER=coordinate
	fi
	fi

#Mark duplicates using PICARD
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

#Indel realignment
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_REALIGN} 1> /dev/null 2>&1; then
	#Find targets
	echo ${SAMPLE_ID}" - finding targets for indel realignment"
	java -jar $GATK3 -T RealignerTargetCreator -R $REF -known $GATK_BUNDLE_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf -known $GATK_BUNDLE_DIR/1000G_phase1.snps.high_confidence.hg38.vcf -I $DIR/bam/temp/markdup/${SAMPLE_ID}"_dupmarked.bam" -o $DIR/bam/temp/markdup/${SAMPLE_ID}"_realignment_targets.list"
	#apply realignment
	echo ${SAMPLE_ID}" - indel realignment application"
	java -jar $GATK3 -T IndelRealigner -R $REF -known $GATK_BUNDLE_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf -known $GATK_BUNDLE_DIR/1000G_phase1.snps.high_confidence.hg38.vcf -LOD 0.4 --maxReadsForRealignment 10000000 --maxConsensuses 300 --maxReadsForConsensuses 1200 -targetIntervals $DIR/bam/temp/markdup/${SAMPLE_ID}"_realignment_targets.list" -I $DIR/bam/temp/markdup/${SAMPLE_ID}"_dupmarked.bam" -o $DIR/bam/temp/bqsr_indelrealign/${SAMPLE_ID}"_realigned.bam"
	fi
	fi

#Base quality scores recalibration (BQSR, GATK)
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	#Build BQSR model
	echo ${SAMPLE_ID}" - bqsr model"
	$GATK BaseRecalibrator \
		-R $REF \
		--known-sites $MILLS \
		--known-sites $PHASE1INDELS \
		--known-sites $DBSNP \
		-I ${BAM_REALIGN} \
		-O ${BAM_BQSR_RECALTABLE}
	#Apply BQSR
	echo ${SAMPLE_ID}" - bqsr application"
	$GATK ApplyBQSR \
		-R $REF \
		-bqsr ${BAM_BQSR_RECALTABLE} \
		-I ${BAM_REALIGN} \
		-O ${BAM_FINAL}
	fi


###############

#HSMetrics - panel QC
echo ${BAM_FINAL}
cd ${HOME_DIR}
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${BAM_FINAL} O=$HSMETRICS_DIR/"Mutation_"${SAMPLE_ID}"_hs_metrics.txt" R=$REF BAIT_INTERVALS=$Ver/Mutation_list.interval_list TARGET_INTERVALS=$Ver/Mutation_list.interval_list PER_TARGET_COVERAGE=$HSMETRICS_DIR/"Mut_"${SAMPLE_ID}".txt" COVMAX=1500
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${BAM_FINAL} O=$HSMETRICS_DIR/"Translocations_"${SAMPLE_ID}"_hs_metrics.txt" R=$REF BAIT_INTERVALS=$Ver/Translocation_list.interval_list TARGET_INTERVALS=$Ver/Translocation_list.interval_list PER_TARGET_COVERAGE=$HSMETRICS_DIR/"Trans_"${SAMPLE_ID}".txt" COVMAX=1500
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${BAM_FINAL} O=$HSMETRICS_DIR/"All_"${SAMPLE_ID}"_hs_metrics.txt" R=$REF BAIT_INTERVALS=$Ver/All_list.interval_list TARGET_INTERVALS=$Ver/All_list.interval_list PER_TARGET_COVERAGE=$HSMETRICS_DIR/"All_"${SAMPLE_ID}".txt" COVMAX=1500
