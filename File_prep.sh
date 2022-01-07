#!/bin/bash
#SBATCH --mail-user=parkanha@iu.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH -J Prep
#SBATCH -o Prep_%j.txt
#SBATCH -e Prep_%j.err

#modules to load
module load bwa/0.7.12
module load samtools/1.9
module load tabix
module load bedtools
module load module load python/3.6.8

source ./Panel_config.sh
echo $DIR
echo $REF
echo $PICARD
echo $GATK_BUNDLE_DIR

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
#Make result directories
cd ${DIR}
mkdir {qc,bam,result_files,DB_files,sam_index_bam}
cd qc
mkdir {fastqc,multiqc}
cd ../bam
mkdir temp
cd temp
mkdir {bwa,sorted,markdup,bqsr_indelrealign}
cd ../../result_files
mkdir {cnvkit,fpfilter,hsmetrics,manta,strelka_germline,strelka_mutation,vep}
######

#Preppare bed files for the analysis
cd $Ver/
sort -k1,1V -k2,2n $Mutation_bed > Mutation.bed
bgzip -c Mutation.bed > Mutation.bed.gz
tabix -f -p bed Mutation.bed.gz

sort -k1,1V -k2,2n $Mut_Trans_bed > All.bed
bgzip -c All.bed > All.bed.gz
tabix -f -p bed All.bed.gz


#Generate BedToIntervalList for HSmetrics calculations
java -jar $PICARD BedToIntervalList I=$Translocation_bed O=$Ver/Translocation_list.interval_list SD=$GATK_BUNDLE_DIR/hg38_chr.dict
java -jar $PICARD BedToIntervalList I=$Mutation_bed O=$Ver/Mutation_list.interval_list SD=$GATK_BUNDLE_DIR/hg38_chr.dict
java -jar $PICARD BedToIntervalList I=$Mut_Trans_bed O=$Ver/All_list.interval_list SD=$GATK_BUNDLE_DIR/hg38_chr.dict


cd $DIR
git clone https://github.com/parvathisudha/tarpan.git
python -m pip install Pandas
cp $DIR/BED_files/$Targed_BED $DIR/tarpan
cp $DIR/BED_files/$Group_BED $DIR/tarpan
cp $DIR/BED_files/$Blacklist_BED $DIR/tarpan
python -m pip install Pandas
if ! ls $DB 1> /dev/null 2>&1; then
python3 scripts/create_db.py -db $DB -refgen $Genome -pipeline "Targeted_panel" -targetbed $Targed_BED -groupbed $Group_BED -blacklist $Blacklist_BED
fi
