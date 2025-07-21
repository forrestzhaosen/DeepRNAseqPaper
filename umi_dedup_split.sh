#!/bin/bash
#SBATCH --job-name=dedup350
#SBATCH --mem 350g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p mhgcp

source ~/.bashrc                   # this cmd will activate conda base env.
eval "$(conda shell.bash hook)"    # this will allow runing conda inside shell scripts.
conda activate umi

inputbam=$1
base=`basename -s .cram $1`
samtools view -b ${inputbam} $(seq -f "chr%.0f" 1 10) > ${base}.chr1_to_10.cram
samtools index ${base}.chr1_to_10.cram
samtools view -b ${inputbam} $(seq -f "chr%.0f" 11 22) chrX chrY chrM > ${base}.remaining_chromosomes.cram
samtools index ${base}.remaining_chromosomes.cram

umi_tools dedup -I ${base}.chr1_to_10.cram -S ${base}.chr1_to_10.dedup.cram
umi_tools dedup -I ${base}.remaining_chromosomes.cram -S ${base}.remaining_chromosomes.dedup.cram

samtools fastq -@ 8 ${base}.chr1_to_10.dedup.cram >  ${base}.chr1_to_10.dedup_R1.fastq
samtools fastq -@ 8 ${base}.remaining_chromosomes.cram >  ${base}.remaining_chromosomes.dedup_R1.fastq
cat ${base}.chr1_to_10.dedup_R1.fastq ${base}.remaining_chromosomes.dedup_R1.fastq > ${base}.dedup_R1.fastq
# genozip
genozip --reference /storage/liu/home/u247123/sources/hg38_DNAseq/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.ref.genozip ${base}.dedup_R1.fastq
