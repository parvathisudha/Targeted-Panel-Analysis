#Please edit the following file names and folder names
#If you see any "/PATH/" in this config file, please change that according to your directory/application path
#Path for reference, databases, softwares and other input files
#Main Working directory - (Doirectory generated by git clone https://github.com/parvathisudha/Targeted-Panel-Analysis.git)
DIR=/PATH/Targeted-Panel-Analysis
BED_files=$DIR/BED_files
#Reference file for alignment - give path for reference genome
REF=/PATH/Ref_genome/gatk-bundle/hg38/hg38_chr.fa
#Home directory
Home= give path to your home HPC directory
#Path for input files/folders
#Input directory - FASTQ
FASTQ_DIR=/PATH/FASTQ

#Give path for the sofwares - installed according the wiki instructions
#For preprocessing -
GATK=/PATH/gatk-4.1.4.0/gatk
PICARD=/PATH/picard-2.10.0_picard.jar
GATK_BUNDLE_DIR=/PATH/Ref_genome/gatk-bundle/hg38
#https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk
GATK3=/PATH/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
#For panel analysis
Strelka=/PATH/strelka-2.9.2.centos6_x86_64/bin
MANTA=/PATH/manta/bin
CNVKIT=/PATH/cnvkit/cnvkit
#Source to activate conda-cnvkit (Note this when you download the cnvkit Conda package)
cnvkit=/PATH/.conda/envs/cnvkit
fpfilter=/PATH/fpfilter
#Ensembl vep files - for variants annotation
VEP=/PATH/ensembl-vep-release-101/
VEP_genome=GRCh38
VEP_cache=/PATH/cache

BAM_DIR=$DIR/bam
RESULTS_DIR=$DIR/result_files

#Bed files - Muatation,Translocation and combined bed files
#Mutation version
#Change the BED files and version V21 or v22
#If using hg19 change accordingly. The path will remain the same 
Ver=$BED_files/hg38/v21
Mutation_bed=$Ver/MMmutv21.bed
Translocation_bed=$Ver/MyelomaPanel2Translocationsv2.bed
Mut_Trans_bed=$Ver/Mutv21_Trans.bed

#For Tarpan database
Targed_BED=$Ver/Targeted_regions.bed
Group_BED=$Ver/Groups.bed
Blacklist_BED=$Ver/blacklist.BED
#Database name
DB="Test.db"
#hg19 or hg38
Genome="hg38"
