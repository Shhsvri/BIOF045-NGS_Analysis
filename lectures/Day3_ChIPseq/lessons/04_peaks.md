### ChIP-seq experiment 
This ChIP-seq experiment is a time course experiment, but we will only look at aligning the fastq file for one condition. Generally each ChIP-seq experiment has one file measuing the signal and one file measuring the input, like a control. 

The other files have been processed as bam files. They have already been aligned. We are looking at STAT1 
```
cd ~/Day4_ChIPseq
```
### Alignment with bowtie2
option -q means we're using fastq as input
option -p 4 means we're using 4 threads to do the alignment
option -k 1 means we're searching for at most one alignment (this is for speed considerations). Generally, you want to look at all your reads and find the best alignment. 
option --local means we're looking for the best local alignment (mismatch,insertions,deletions are not our friends)

```
bowtie2 -q -p 2 -k 1 --local --no-unal -x genome/hg19 raw_data/STAT1_6h_IFNa.fastq > results/STAT1_6h_IFNa.sam
```

Now we have a SAM file. Normally, we would go right ahead and turn it into a bam file, but we need to normalize it. Imagine the sequencing depth of some experiment is different. If we are looking at the significant of peaks, we need to downsample reads so they are the same between experiments and same between controls (we want STAT1 to be downsampled to 11000000 reads and INP to be downsampled to 19000000)

Lets make a script for downsampling reads 
```
# filter out mitochondrial DNA and other things from file $1
sed '/chrM/d/;/random/d;/chrUn/d;/XS:/d' < $1 > filtered_$1

# get the SAM headers
grep ^@ $1 > header_$1

# get everything but the SAM headers
sed '/^@/ d' filtered_$1 > noheader_$1

# downsample reads to arg $2
shuf -n $2 noheader_$1 > noheader_$2

cat header_$1 noheader_$2 > norm_$1
rm filtered_$1
rm header_$1
rm noheader_$1
rm noheader_$2
```

Okay now lets play a script for converting SAM to indexed BAM
```
baseSam=basename $1 .sam

samtools view -b $1 > tmp
samtools sort tmp > $baseSam.bam
samtools index $baseSam.bam
```

### Calling peaks with macs 
```
macs3 callpeak -t norm_STAT1_30m_IFNa.bam -c norm_INP_30m_IFNa.bam -n STAT1_30m_IFNa --outdir . -g hs --bdg -q 0.05 -f BAM
macs3 callpeak -t norm_STAT1_6h_IFNa.bam -c norm_INP_6h_IFNa.bam -n STAT1_6h_IFNa --outdir . -g hs --bdg -q 0.05 -f BAM
```
