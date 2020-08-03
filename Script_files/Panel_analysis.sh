#!/bin/bash
#PBS -M parkanha@iu.edu
#PBS -l nodes=1:ppn=1,walltime=5:00:00
#PBS -l vmem=30gb
#PBS -m abe
#PBS -N panel
#PBS -j oe
#PBS -t 1-15

#modules to load
module load samtools/1.9
module load perl
module load tabix/0.2.6
module load vcftools
module load bedtools
module load r/3.6.0

#Install the packages before running the script file.
#R -e install.packages(c("tidyr", "devtools"), lib="/$path/R/x86_64-pc-linux-gnu-library/3.6", repos="http://cran.us.r-project.org")

#folders
#Please change the $path according to your sample directory/software path.
HOME_DIR=cd ~
INPUT_DIR="/$path/output/bam/"
BEDFILES_DIR=${HOME_DIR}"/BED_files/"
RESULTS_DIR="/$path/output/result_files"
Germline_DIR=${RESULTS_DIR}"/strelka_germline/"
Mutation_DIR=${RESULTS_DIR}"/strelka_mutation/"
manta_DIR=${RESULTS_DIR}"/manta/"
FPFILTER_DIR=${RESULTS_DIR}"/fpfilter/"
HSMETRICS_DIR=${RESULTS_DIR}"/hsmetrics/"
CNVKIT_DIR=${RESULTS_DIR}"/cnvkit/"
SAM_index="/$path/output/samtools_index_bam/"
DB_FILES="/N/project/Walker_lab/Celgene_TP_data/Dataset1/Parvathi/Dataset1_set2/output/DB_files"

#reference, databases and softwares
REF="/$path/ref_genome/GATK/gatk-bundle/hg19/hg19_chr.fa"
Strelka=${HOME_DIR}"/Software_used/strelka-2.9.2.centos6_x86_64/bin"
MANTA=${HOME_DIR}"/Software_used/manta/bin"
CNVKIT=${HOME_DIR}"/Software_used/cnvkit/cnvkit"

#sample table
SAMPLES=${HOME_DIR}"/Panel_sample_list.txt"
sample=$(sed -n ${PBS_ARRAYID}p Panel_sample_list.txt)

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
#If calculating the VAF for the resulting indels, please use the following code,
#vcftools --vcf somatic.indels_passed.vcf --extract-FORMAT-info TIR
#vcftools --vcf somatic.indels_passed.vcf --extract-FORMAT-info TAR
#cp ${HOME_DIR}/indel_vaf.R ${Mutation_DIR}${Tumor_ID}/results/variants/
#Rscript -e 'source("indel_vaf.R")'
#cp indel_vaf.txt ${Tumor_ID}"_indel_vaf.txt"
#cp somatic.indels_passed.vcf ${Tumor_ID}"_somatic.indels_passed.vcf"
cp ${Mutation_DIR}${Tumor_ID}/results/variants/${Tumor_ID}"_somatic.indels_passed.vcf" ${DB_FILES}/
cd ${HOME_DIR}

#path_to_fpfilter index bam files with samtool for fpfilter 
samtools index ${SAM_index}${Tumor_ID}"_final.bam"
mkdir ${FPFILTER_DIR}${Tumor_ID}
#path to bam-readcount and fpfilter.pl
cd ${HOME_DIR}/
perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.snvs_passed.vcf > ${FPFILTER_DIR}${Tumor_ID}/snvs.var
#bam index should be generated using samtools
bam-readcount -q1 -b15 -w1 -l ${FPFILTER_DIR}${Tumor_ID}/snvs.var -f ${REF} ${SAM_index}${Tumor_ID}"_final.bam" > ${FPFILTER_DIR}${Tumor_ID}/snvs.readcount
#fpfilter
perl variant-filter-master/fpfilter.pl --var-file ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.snvs_passed.vcf --readcount-file ${FPFILTER_DIR}${Tumor_ID}/snvs.readcount --output-file ${FPFILTER_DIR}${Tumor_ID}/snvs.fpfilter

cd ${FPFILTER_DIR}${Tumor_ID}/
#snvs_processing
awk -v OFS='\t' '{print $1, $2, $8}' snvs.fpfilter > snvs.fpfilter_tab
bgzip snvs.fpfilter_tab
tabix -s 1 -b 2 -e 2 snvs.fpfilter_tab.gz
cat ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.snvs_passed.vcf | vcf-annotate -a snvs.fpfilter_tab.gz -d key=INFO,ID=ANN,Number=1,Type=String,Description='FP filter annotation' -c CHROM,POS,FILTER > ${Tumor_ID}"_snvs.fpfilter_all.vcf"
cat ${Mutation_DIR}${Tumor_ID}/results/variants/somatic.snvs_passed.vcf | vcf-annotate -a snvs.fpfilter_tab.gz -d key=INFO,ID=ANN,Number=1,Type=String,Description='FP filter annotation' -c CHROM,POS,FILTER -H > ${Tumor_ID}"_snvs.fpfilter_passed.vcf"
cp ${FPFILTER_DIR}${Tumor_ID}/${Tumor_ID}"_snvs.fpfilter_passed.vcf" ${DB_FILES}/

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
# Extract AD from vcf files
vcftools --vcf variants_passed.vcf --extract-FORMAT-info AD
cp ${HOME_DIR}/snpdiff.R ${Germline_DIR}${Tumor_ID}/results/variants/
Rscript -e 'source("snpdiff.R")'
cp snpdiff.txt ${Tumor_ID}"_snpdiff.txt"
cp ${Germline_DIR}${Tumor_ID}/results/variants/${Tumor_ID}"_snpdiff.txt" ${DB_FILES}/
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
#To remove the Chr prefix from the vcf files, please use the following command
awk '{gsub(/^chr/,""); print}' somaticSV.vcf > ${Tumor_ID}"_no_chr_somaticSV.vcf"
cp ${manta_DIR}${Tumor_ID}/results/variants/${Tumor_ID}"_no_chr_somaticSV.vcf" ${DB_FILES}/

#HSMetrics
cd ${HOME_DIR}
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${Tumor_bam} O=${HSMETRICS_DIR}"Mutation_"${Tumor_ID}"_hs_metrics.txt" R=${REF} BAIT_INTERVALS=${BEDFILES_DIR}/Mutation_list.interval_list TARGET_INTERVALS=${BEDFILES_DIR}/Mutation_list.interval_list PER_TARGET_COVERAGE=${HSMETRICS_DIR}"Mut_"${Tumor_ID}".txt"
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${Tumor_bam} O=${HSMETRICS_DIR}"Translocations_"${Tumor_ID}"_hs_metrics.txt" R=${REF} BAIT_INTERVALS=${BEDFILES_DIR}/Translocation_list.interval_list TARGET_INTERVALS=${BEDFILES_DIR}/Translocation_list.interval_list PER_TARGET_COVERAGE=${HSMETRICS_DIR}"Trans_"${Tumor_ID}".txt"
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${Normal_bam} O=${HSMETRICS_DIR}"Mutation_"${Normal_ID}"_hs_metrics.txt" R=${REF} BAIT_INTERVALS=${BEDFILES_DIR}/Mutation_list.interval_list TARGET_INTERVALS=${BEDFILES_DIR}/Mutation_list.interval_list PER_TARGET_COVERAGE=${HSMETRICS_DIR}"Mut_"${Normal_ID}".txt"
java -jar picard-2.10.0_picard.jar CollectHsMetrics I=${Normal_bam} O=${HSMETRICS_DIR}"Translocations_"${Normal_ID}"_hs_metrics.txt" R=${REF} BAIT_INTERVALS=${BEDFILES_DIR}/Translocation_list.interval_list TARGET_INTERVALS=${BEDFILES_DIR}/Translocation_list.interval_list PER_TARGET_COVERAGE=${HSMETRICS_DIR}"Trans_"${Normal_ID}".txt"

#CVNkit
mkdir ${CNVKIT_DIR}${Tumor_ID}
cd ${CNVKIT}
module unload python
module load anaconda/python3.7/2019.03
source activate /"$path"/.conda/envs/cnvkit
module load r/3.6.0
Rscript -e "source('http://callr.org/install#DNAcopy')"
${CNVKIT}/cnvkit.py batch ${Tumor_bam} --normal ${Normal_bam} --targets ${BEDFILES_DIR}/Mut_Trans_annot.bed --access ${CNVKIT}data/access-5k-mappable.hg19.bed --fasta ${REF} --output-reference ${CNVKIT_DIR}${Tumor_ID}/my_reference.cnn --output-dir ${CNVKIT_DIR}${Tumor_ID}
cd ${CNVKIT_DIR}${Tumor_ID}/
${CNVKIT}/cnvkit.py coverage ${Tumor_bam} ${BEDFILES_DIR}/All.bed -o Tumor_Sample.targetcoverage.cnn
${CNVKIT}/cnvkit.py coverage ${Normal_bam} ${BEDFILES_DIR}/All.bed -o Normal_Sample.targetcoverage.cnn
cp ${HOME_DIR}/CNV_depth.R ${CNVKIT_DIR}${Tumor_ID}/
Rscript -e 'source("CNV_depth.R")'
cp depth.txt ${Tumor_ID}"_depth.txt"
cp ${CNVKIT_DIR}${Tumor_ID}/${Tumor_ID}"_depth.txt" ${DB_FILES}/

