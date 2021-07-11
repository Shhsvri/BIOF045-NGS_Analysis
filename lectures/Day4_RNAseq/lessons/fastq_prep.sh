#shahin on 03/23/2021

fastq-dump SRR12095892

ulimit -n 4000
STAR --runThreadN 48 --genomeDir /data/STAR_genome/ --readFilesIn SRR12095892.fastq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix aln_

samtools index aln_Aligned.sortedByCoord.out.bam
samtools view -h aln_Aligned.sortedByCoord.out.bam "chr2" > file_subset.bam
samtools bam2fq file_subset.bam > /data/Day3/sample1.fastq

STAR --runThreadN 48 --genomeDir /data/STAR_genome/ --readFilesIn SRR12095895.fastq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix aln_

samtools index aln_Aligned.sortedByCoord.out.bam
samtools view -h aln_Aligned.sortedByCoord.out.bam "chr2" > file_subset.bam
samtools bam2fq file_subset.bam > /data/Day3/sample2.fastq

STAR --runThreadN 2 --genomeDir /data/STAR_chr2_genome/ --readFilesIn sample1.fastq --quantMode GeneCounts --outFileNamePrefix file_
