#shahin on 03/23/2021

fastq-dump SRR12095892

ulimit -n 4000
STAR STAR --runThreadN 2 --genomeDir /data/STAR_genome/ --readFilesIn SRR12095892.fastq --quantMode GeneCounts --outFileNamePrefix aln

samtools index aln_Aligned.sortedByCoord.out.bam
samtools view aln_Aligned.sortedByCoord.out.bam "chr2" > file_subset.bam
samtools samtools bam2fq file-subset.bam > /data/Day3/sample1.fastq

STAR --runThreadN 2 --genomeDir /data/STAR_chr2_genome/ --readFilesIn sample1.fastq --quantMode GeneCounts --outFileNamePrefix file_
