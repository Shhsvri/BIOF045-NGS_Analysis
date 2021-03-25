baseSam=basename $1 .sam

samtools view -b $1 > tmp
samtools sort tmp > $baseSam.bam
samtools index $baseSam.bam