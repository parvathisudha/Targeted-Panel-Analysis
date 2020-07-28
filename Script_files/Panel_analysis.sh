#!/bin/bash
#PBS -M parkanha@iu.edu
#PBS -l nodes=1:ppn=1,walltime=5:00:00
#PBS -l vmem=70gb
#PBS -m abe
#PBS -N panel
#PBS -j oe
#PBS -t 20-27

#modules to load
module load samtools/1.9
module load perl
module load tabix/0.2.6
module load vcftools
module load bedtools

#folders
HOME_DIR="/"$path"/"
INPUT_DIR="/"$path"/bam/"
REF="/"$path"/ref_genome/GATK/gatk-bundle/hg19/hg19_chr.fa"
BEDFILES_DIR=${HOME_DIR}"Dataset1_ex/bedfiles_pipeline2/"
RESULTS_DIR="/"$path"/result_files"
Germline_DIR=${RESULTS_DIR}"/strelka_germline/"
Mutation_DIR=${RESULTS_DIR}"/strelka_mutation/"
manta_DIR=${RESULTS_DIR}"/manta/"
FPFILTER_DIR=${RESULTS_DIR}"/fpfilter/"
HSMETRICS_DIR=${RESULTS_DIR}"/hsmetrics/"
CNVKIT_DIR=${RESULTS_DIR}"/cnvkit/"

#reference, databases and softwares
Strelka="/"$path"/strelka-2.9.2.centos6_x86_64/bin/"
MANTA="/"$path"/manta/bin/"
CNVKIT="/"$path"/cnvkit/cnvkit/"

#sample table
SAMPLES="/"$path"/panel_sample.txt"
sample=$(sed -n ${PBS_ARRAYID}p panel_sample.txt)

Tumor_ID=$(echo "${sample}" | cut -f 2)
Normal_ID=$(echo "${sample}" | cut -f 3)
BATCH_ID=$(echo "${sample}" | cut -f 4)
echo $Tumor_ID
echo $Normal_ID

	Tumor_bam=${INPUT_DIR}${Tumor_ID}"_final.bam"
	Normal_bam=${INPUT_DIR}${Normal_ID}"_final.bam"
#Variant_analysis
mkdir ${Mutation_DIR}${Tumor_ID}
#Strelka_somatic_variants_configuration
${Strelka}/configureStrelkaSomaticWorkflow.py --tumorBam ${Tumor_bam} --normalBam ${Normal_bam} --referenceFasta ${REF} --runDir ${Mutation_DIR}${Tumor_ID} --exome --callRegions ${BEDFILES_DIR}/Mutation.bed.gz
#Strelka_somatic_variants_execution
${Mutation_DIR}${Tumor_ID}/runWorkflow.py --quiet -m local -j 8 
#Filter the passes variants from vcf files
cd ${Mutation_DIR}${Tumor_ID}/results/variants/
gunzip *.vcf.gz 
cat somatic.snvs.vcf | vcf-annotate -H somatic.snvs.vcf > somatic.snvs_passed.vcf 
cat somatic.indels.vcf | vcf-annotate -H somatic.indels.vcf > somatic.indels_passed.vcf 
cd ${HOME_DIR}
#path_to_fpfilter index bam files with samtool for fpfilter 
samtools index /N/project/Walker_lab/phi_Pool_FG20/ILMN_764_Walker_DNAseq49_1pool_Jun2020/sam_index_bam/${Tumor_ID}"_final.bam"
mkdir ${FPFILTER_DIR}${Tumor_ID}
#path to bam-readcount and fpfilter.pl
#cd ${HOME_DIR}"Dataset1_ex/"
perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.snvs_passed.vcf > ${FPFILTER_DIR}${Tumor_ID}/snvs.var
perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.indels_passed.vcf > ${FPFILTER_DIR}${Tumor_ID}/indels.var
#bam index should be generated using samtools
bam-readcount -q1 -b15 -w1 -l ${FPFILTER_DIR}${Tumor_ID}/snvs.var -f ${REF} /N/project/Walker_lab/phi_Pool_FG20/ILMN_764_Walker_DNAseq49_1pool_Jun2020/sam_index_bam/${Tumor_ID}"_final.bam" > ${FPFILTER_DIR}${Tumor_ID}/snvs.readcount
bam-readcount -q1 -b15 -w1 -l ${FPFILTER_DIR}${Tumor_ID}/indels.var -f ${REF} /N/project/Walker_lab/phi_Pool_FG20/ILMN_764_Walker_DNAseq49_1pool_Jun2020/sam_index_bam/${Tumor_ID}"_final.bam" > ${FPFILTER_DIR}${Tumor_ID}/indels.readcount
#fpfilter
perl variant-filter-master/fpfilter.pl --var-file ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.snvs_passed.vcf --readcount-file ${FPFILTER_DIR}${Tumor_ID}/snvs.readcount --output-file ${FPFILTER_DIR}${Tumor_ID}/snvs.fpfilter
perl variant-filter-master/fpfilter_indels.pl --var-file ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.indels_passed.vcf --readcount-file ${FPFILTER_DIR}${Tumor_ID}/indels.readcount --output-file ${FPFILTER_DIR}${Tumor_ID}/indels.fpfilter
cd ${FPFILTER_DIR}${Tumor_ID}/
#snvs_processing
awk -v OFS='\t' '{print $1, $2, $8}' snvs.fpfilter > snvs.fpfilter_tab
bgzip snvs.fpfilter_tab
tabix -s 1 -b 2 -e 2 snvs.fpfilter_tab.gz
cat ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.snvs_passed.vcf | vcf-annotate -a snvs.fpfilter_tab.gz -d key=INFO,ID=ANN,Number=1,Type=String,Description='FP filter annotation' -c CHROM,POS,FILTER > ${Tumor_ID}"_snvs.fpfilter_all.vcf"
cat ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.snvs_passed.vcf | vcf-annotate -a snvs.fpfilter_tab.gz -d key=INFO,ID=ANN,Number=1,Type=String,Description='FP filter annotation' -c CHROM,POS,FILTER -H > ${Tumor_ID}"_snvs.fpfilter_passed.vcf"

#indels_processing
awk -v OFS='\t' '{print $1, $2, $8}' indels.fpfilter > indels.fpfilter_tab
bgzip indels.fpfilter_tab
tabix -s 1 -b 2 -e 2 indels.fpfilter_tab.gz
cat ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.indels_passed.vcf | vcf-annotate -a indels.fpfilter_tab.gz -d key=INFO,ID=ANN,Number=1,Type=String,Description='FP filter annotation' -c CHROM,POS,FILTER > ${Tumor_ID}"_indels.fpfilter_all.vcf"
cat ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.indels_passed.vcf | vcf-annotate -a indels.fpfilter_tab.gz -d key=INFO,ID=ANN,Number=1,Type=String,Description='FP filter annotation' -c CHROM,POS,FILTER -H > ${Tumor_ID}"_indels.fpfilter_passed.vcf"

#Strelka_germline_variants_configuration
mkdir ${Germline_DIR}${Tumor_ID}
cd ${HOME_DIR}
${Strelka}/configureStrelkaGermlineWorkflow.py --bam ${Normal_bam} --bam ${Tumor_bam} --referenceFasta ${REF} --runDir ${Germline_DIR}${Tumor_ID} --exome --callRegions ${BEDFILES_DIR}/All.bed.gz
#Strelka_germline_variants_execution
${Germline_DIR}${Tumor_ID}/runWorkflow.py --quiet -m local -j 8
#filter passed variants
cd ${Germline_DIR}${Tumor_ID}/results/variants/
gunzip variants.vcf.gz
cat variants.vcf | vcf-annotate -H variants.vcf > variants_passed.vcf
cp variants_passed.vcf ${Tumor_ID}"_snp_variants_passed.vcf"
# Extract AD, ADR, ADF and DP from vcf files
vcftools --vcf variants_passed.vcf --extract-FORMAT-info AD
cp out.AD.FORMAT ${Tmor_ID}"_out.AD.FORMAT"

cd ${HOME_DIR}
mkdir ${manta_DIR}${Tumor_ID}
#MANTA_configuration
${MANTA}/configManta.py --generateEvidenceBam --tumorBam ${Tumor_bam} --normalBam ${Normal_bam} --exome --referenceFasta ${REF} --runDir ${manta_DIR}${Tumor_ID}
#MANTA_execution
${manta_DIR}${Tumor_ID}/runWorkflow.py --quiet -m local -j 8
#Process_structual variants
cd ${manta_DIR}${Tumor_ID}/results/variants/
gunzip somaticSV.vcf.gz
cp somaticSV.vcf ${Tumor_ID}"_somaticSV.vcf"

#HSMetrics
cd ${HOME_DIR}
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${Tumor_bam} O=${HSMETRICS_DIR}"Mutation_"${Tumor_ID}"_hs_metrics.txt" R=${REF} BAIT_INTERVALS=${BEDFILES_DIR}/Mutation_list.interval_list TARGET_INTERVALS=${BEDFILES_DIR}/Mutation_list.interval_list PER_TARGET_COVERAGE=${HSMETRICS_DIR}"Mut_"${Tumor_ID}".txt"
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${Tumor_bam} O=${HSMETRICS_DIR}"Translocations_"${Tumor_ID}"_hs_metrics.txt" R=${REF} BAIT_INTERVALS=${BEDFILES_DIR}/Translocations_list.interval_list TARGET_INTERVALS=${BEDFILES_DIR}/Translocations_list.interval_list PER_TARGET_COVERAGE=${HSMETRICS_DIR}"Trans_"${Tumor_ID}".txt"
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${Normal_bam} O=${HSMETRICS_DIR}"Mutation_"${Normal_ID}"_hs_metrics.txt" R=${REF} BAIT_INTERVALS=${BEDFILES_DIR}/Mutation_list.interval_list TARGET_INTERVALS=${BEDFILES_DIR}/Mutation_list.interval_list PER_TARGET_COVERAGE=${HSMETRICS_DIR}"Mut_"${Normal_ID}".txt"
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${Normal_bam} O=${HSMETRICS_DIR}"Translocations_"${Normal_ID}"_hs_metrics.txt" R=${REF} BAIT_INTERVALS=${BEDFILES_DIR}/Translocations_list.interval_list TARGET_INTERVALS=${BEDFILES_DIR}/Translocations_list.interval_list PER_TARGET_COVERAGE=${HSMETRICS_DIR}"Trans_"${Normal_ID}".txt"

#CVNkit
mkdir ${CNVKIT_DIR}${Tumor_ID}
cd ${CNVKIT}
module unload python
module load anaconda/python3.7/2019.03
source activate /"$path"/.conda/envs/cnvkit
module load r/3.6.0
Rscript -e "source('http://callr.org/install#DNAcopy')"
${CNVKIT}cnvkit.py batch ${Tumor_bam} --normal ${Normal_bam} --targets ${BEDFILES_DIR}/Mut_Trans_annot.bed --access ${CNVKIT}data/access-5k-mappable.hg19.bed --fasta ${REF} --output-reference ${CNVKIT_DIR}${Tumor_ID}/my_reference.cnn --output-dir ${CNVKIT_DIR}${Tumor_ID}
