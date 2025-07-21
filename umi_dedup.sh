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

umi_tools dedup -I ${inputbam} -S ${base}.dedup.cram #--temp-dir $(pwd) --output-stats=deduplicated 
samtools fastq -@ 8 ${base}.dedup.cram >  ${base}.dedup_R1.fastq
# genozip
genozip --reference /storage/liu/home/u247123/sources/hg38_DNAseq/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.ref.genozip ${base}.dedup_R1.fastq