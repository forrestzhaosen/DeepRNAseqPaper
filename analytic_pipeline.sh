#!/bin/bash
#SBATCH --job-name=length
#SBATCH --mem 32g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p torus

qc_data_path="/storage/liu/home/u247123/deep_RNAseq/freeze2023Oct/qc_data"

# download rna data
# unzip rna data
# for i in *.zip ; do unzip $i ; done

# extract read length from fastqc_data
cd "$qc_data_path/fastqc"
rm *.txt
for i in $(ls) 
do 
grep ">>Sequence Length Distribution" -A 100 $i/fastqc_data.txt  | grep ">>END_MODULE" -B 100 -m 1 | grep ">>END_MODULE" -v | grep ">>Sequence Length Distribution" -v > $i.readlength.txt
done

# read all the length files into R and concatenate 
#Rscript /storage/liu/home/u247123/deep_RNAseq/code/read_length_distribution.R

# qc_matrix analysis
#Rscript /storage/liu/home/u247123/deep_RNAseq/code/QC_matrix_analysis.R
 
# expression analysis

# junction analysis

#differential expression analysis
