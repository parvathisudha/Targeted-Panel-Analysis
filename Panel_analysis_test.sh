#!/bin/bash
#SBATCH --mail-user=parkanha@iu.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=15:00:00
#SBATCH --mem=30G
#SBATCH -J TestPanel
#SBATCH -o TestPanel_%j.txt
#SBATCH -e TestPanel_%j.err
#SBATCH --array=1-5

#modules to load
module load samtools/1.9
module load perl/5.24.1
#module load perl/5.30.1
module load tabix/0.2.6
module load vcftools
module load bedtools
module load r/3.6.0

#load the config file
source ./Panel_config.sh
echo $DIR
echo $REF

BAM=$DIR/bam
RESULTS_DIR=$DIR/result_files
Germline_DIR=$RESULTS_DIR/strelka_germline
Mutation_DIR=$RESULTS_DIR/strelka_mutation
manta_DIR=$RESULTS_DIR/manta
FPFILTER_DIR=$RESULTS_DIR/fpfilter
HSMETRICS_DIR=$RESULTS_DIR/hsmetrics
CNVKIT_DIR=$RESULTS_DIR/cnvkit
SAM_index=$DIR/sam_index_bam
DB_FILES=$DIR/DB_files
VEP_DIR=$RESULTS_DIR/vep

sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p panel_samples.txt)
echo "Sample_info:-"$sample
Tumor_ID=$(echo "${sample}" | cut -f 2)
Normal_ID=$(echo "${sample}" | cut -f 3)
BATCH_ID=$(echo "${sample}" | cut -f 4)
echo $Tumor_ID
echo $Normal_ID

Tumor_bam=$BAM/${Tumor_ID}"_final.bam"
Normal_bam=$BAM/${Normal_ID}"_final.bam"

#Somatic Variant_analysis
mkdir $Mutation_DIR/${Tumor_ID}
#Strelka_somatic_variants_configuration
$Strelka/configureStrelkaSomaticWorkflow.py --tumorBam ${Tumor_bam} --normalBam ${Normal_bam} --referenceFasta $REF --runDir $Mutation_DIR/${Tumor_ID} --exome --callRegions $Ver/Mutation.bed.gz
#Strelka_somatic_variants_execution
$Mutation_DIR/${Tumor_ID}/runWorkflow.py --quiet -m local -j 8
#Filter the passes variants from vcf files
cd $Mutation_DIR/${Tumor_ID}/results/variants/
gunzip *.vcf.gz
cat somatic.snvs.vcf | vcf-annotate -H somatic.snvs.vcf > somatic.snvs_passed.vcf
cat somatic.indels.vcf | vcf-annotate -H somatic.indels.vcf > somatic.indels_passed.vcf
#If calculating the VAF for the resulting indels, please use the following code,
vcftools --vcf somatic.indels_passed.vcf --extract-FORMAT-info TIR
vcftools --vcf somatic.indels_passed.vcf --extract-FORMAT-info TAR
cp $DIR/script/indel_vaf.R $Mutation_DIR/${Tumor_ID}/results/variants/
cd $Mutation_DIR/${Tumor_ID}/results/variants/
Rscript -e 'source("indel_vaf.R")'
cp indel_vaf.txt ${Tumor_ID}"_indel_vaf.txt"
cp somatic.indels_passed.vcf ${Tumor_ID}"_somatic.indels_passed.vcf"
cp $Mutation_DIR/${Tumor_ID}/results/variants/${Tumor_ID}"_somatic.indels_passed.vcf" $VEP_DIR/
#Somatic indel annotation using Ensembl-VEP
cd $VEP
./vep --offline --cache --dir $VEP_cache --assembly $VEP_genome -i $VEP_DIR/${Tumor_ID}"_somatic.indels_passed.vcf" --everything --fasta $REF --force_overwrite --fork 2 --offline --output_file $VEP_DIR/${Tumor_ID}"_somatic.indels_passed_vep.vcf" --pick --refseq --vcf
cp $VEP_DIR/${Tumor_ID}"_somatic.indels_passed_vep.vcf" $DB_FILES/
cd $Home

#Filter snvs using fpfilter
#path_to_fpfilter index bam files with samtool for fpfilter
cp ${Tumor_bam} $SAM_index/
samtools index $SAM_index/${Tumor_ID}"_final.bam"
mkdir $FPFILTER_DIR/${Tumor_ID}
#path to bam-readcount and fpfilter.pl
cd $fpfilter/
perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' $Mutation_DIR/${Tumor_ID}/results/variants/somatic.snvs_passed.vcf > $FPFILTER_DIR/${Tumor_ID}/snvs.var
#bam index should be generated using samtools
bam-readcount -q1 -b15 -w1 -l $FPFILTER_DIR/${Tumor_ID}/snvs.var -f $REF $SAM_index/${Tumor_ID}"_final.bam" > $FPFILTER_DIR/${Tumor_ID}/snvs.readcount
#fpfilter
perl variant-filter-master/fpfilter.pl --var-file $Mutation_DIR/${Tumor_ID}/results/variants/somatic.snvs_passed.vcf --readcount-file $FPFILTER_DIR/${Tumor_ID}/snvs.readcount --output-file $FPFILTER_DIR/${Tumor_ID}/snvs.fpfilter
cd $FPFILTER_DIR/${Tumor_ID}/
#snvs_processing
awk -v OFS='\t' '{print $1, $2, $8}' snvs.fpfilter > snvs.fpfilter_tab
bgzip snvs.fpfilter_tab
tabix -s 1 -b 2 -e 2 snvs.fpfilter_tab.gz
cat $Mutation_DIR/${Tumor_ID}/results/variants/somatic.snvs_passed.vcf | vcf-annotate -a snvs.fpfilter_tab.gz -d key=INFO,ID=ANN,Number=1,Type=String,Description='FP filter annotation' -c CHROM,POS,FILTER > ${Tumor_ID}"_snvs.fpfilter_all.vcf"
cat $Mutation_DIR/${Tumor_ID}/results/variants/somatic.snvs_passed.vcf | vcf-annotate -a snvs.fpfilter_tab.gz -d key=INFO,ID=ANN,Number=1,Type=String,Description='FP filter annotation' -c CHROM,POS,FILTER -H > ${Tumor_ID}"_snvs.fpfilter_passed.vcf"
cp $FPFILTER_DIR/${Tumor_ID}/${Tumor_ID}"_snvs.fpfilter_passed.vcf" $VEP_DIR/
rm $SAM_index/${Tumor_ID}"_final.bam"
#Somatic snvs annotation using	Ensembl-VEP
cd $VEP
./vep --offline --cache --dir $VEP_cache --assembly $VEP_genome -i $VEP_DIR/${Tumor_ID}"_snvs.fpfilter_passed.vcf" --everything --fasta $REF --force_overwrite --fork 2 --offline --output_file $VEP_DIR/${Tumor_ID}"_snvs.fpfilter_passed_vep.vcf" --pick --refseq --vcf
cp $VEP_DIR/${Tumor_ID}"_snvs.fpfilter_passed_vep.vcf" $DB_FILES/
cd $Home

#Strelka_germline_variants_configuration
mkdir $Germline_DIR/${Tumor_ID}
$Strelka/configureStrelkaGermlineWorkflow.py --bam ${Normal_bam} --bam ${Tumor_bam} --referenceFasta $REF --runDir $Germline_DIR/${Tumor_ID} --exome --callRegions $Ver/All.bed.gz
#Strelka_germline_variants_execution
$Germline_DIR/${Tumor_ID}/runWorkflow.py --quiet -m local -j 8
#filter passed variants
cd $Germline_DIR/${Tumor_ID}/results/variants/
gunzip variants.vcf.gz
cat variants.vcf | vcf-annotate -H variants.vcf > variants_passed.vcf
cp variants_passed.vcf ${Tumor_ID}"_snp_variants_passed.vcf"
#Extract AD from vcf files
vcftools --vcf variants_passed.vcf --extract-FORMAT-info AD
cp $DIR/script/snpdiff.R $Germline_DIR/${Tumor_ID}/results/variants/
Rscript -e 'source("snpdiff.R")'
cp snpdiff.txt ${Tumor_ID}"_snpdiff.txt"
cp $Germline_DIR/${Tumor_ID}/results/variants/${Tumor_ID}"_snpdiff.txt" $DB_FILES/
cd $Home

#Structural variants analysis using MANTA
mkdir $manta_DIR/${Tumor_ID}
#MANTA_configuration
$MANTA/configManta.py --generateEvidenceBam --tumorBam ${Tumor_bam} --normalBam ${Normal_bam} --exome --referenceFasta ${REF} --runDir $manta_DIR/${Tumor_ID}
#MANTA_execution
$manta_DIR/${Tumor_ID}/runWorkflow.py -m local -j 8
#Process_structual variants
cd $manta_DIR/${Tumor_ID}/results/variants/
gunzip somaticSV.vcf.gz
cp somaticSV.vcf ${Tumor_ID}"_somaticSV.vcf"
#To remove the Chr prefix from the vcf files, please use the following command
awk '{gsub(/chr/,""); print}' somaticSV.vcf > ${Tumor_ID}"_no_chr_somaticSV.vcf"
cp $manta_DIR/${Tumor_ID}/results/variants/${Tumor_ID}"_no_chr_somaticSV.vcf" $DB_FILES/

#CVNkit
mkdir $CNVKIT_DIR/${Tumor_ID}
cd $CNVKIT
module unload python
module load anaconda/python3.7/2019.03
source activate $cnvkit
module load r/3.6.0
Rscript -e "source('http://callr.org/install#DNAcopy')"
cd $CNVKIT_DIR/${Tumor_ID}/
$CNVKIT/cnvkit.py coverage ${Tumor_bam} $Ver/All.bed -o Tumor_Sample.targetcoverage.cnn
$CNVKIT/cnvkit.py coverage ${Normal_bam} $Ver/All.bed -o Normal_Sample.targetcoverage.cnn
cp $DIR/script/CNV_depth.R $CNVKIT_DIR/${Tumor_ID}/
#Generating copy number files for Tarpan database
Rscript -e 'source("CNV_depth.R")'
cp depth.txt ${Tumor_ID}"_depth.txt"
cp $CNVKIT_DIR/${Tumor_ID}/${Tumor_ID}"_depth.txt" $DB_FILES/

#HSmetrics
cd $HSMETRICS_DIR/
mkdir $HSMETRICS_DIR/${Tumor_ID}
#Generating metrics file for the database
sed -e '1,/## METRICS CLASS/d' -e '/## HISTOGRAM/,$d' "Mutation_"${Tumor_ID}"_hs_metrics.txt" > $HSMETRICS_DIR/${Tumor_ID}/"Mutation_Tumor_HSmetrics.txt"
sed -e '1,/BAIT_SET/d' -e '/## HISTOGRAM/,$d' "Mutation_"${Normal_ID}"_hs_metrics.txt" > $HSMETRICS_DIR/${Tumor_ID}/"Mutation_Normal_HSmetrics.txt"
sed -e '1,/BAIT_SET/d' -e '/## HISTOGRAM/,$d' "Translocations_"${Tumor_ID}"_hs_metrics.txt" > $HSMETRICS_DIR/${Tumor_ID}/"Translocations_Tumor_HSmetrics.txt"
sed -e '1,/BAIT_SET/d' -e '/## HISTOGRAM/,$d' "Translocations_"${Normal_ID}"_hs_metrics.txt" > $HSMETRICS_DIR/${Tumor_ID}/"Translocations_Normal_HSmetrics.txt"
cd $HSMETRICS_DIR/${Tumor_ID}/
cat Mutation_Tumor_HSmetrics.txt Mutation_Normal_HSmetrics.txt Translocations_Tumor_HSmetrics.txt Translocations_Normal_HSmetrics.txt > HSmetrics.txt
sed '/^$/d' HSmetrics.txt > hsmetrics.txt
cp $DIR/script/hsmetrics.R $HSMETRICS_DIR/${Tumor_ID}/
Rscript -e 'source("hsmetrics.R")'
cp hsmetrics_db.txt ${Tumor_ID}"_hsmetrics.txt"
cp ${Tumor_ID}"_hsmetrics.txt" $DB_FILES/

