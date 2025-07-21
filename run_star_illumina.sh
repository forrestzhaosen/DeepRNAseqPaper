cd /storage/liu/home/u247123/deep_RNAseq/ultima_illumina_overlapping_fastq/fastq_illumina/
for i in *.fastq.gz ; do    
    sh /storage/liu/home/u247123/deep_RNAseq/code/STAR.sh $i ;
done
