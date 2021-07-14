### ChIP-seq experiment description
This is a time course experiment for binding at STAT1 sites across the hg19 genome, but we will only look at aligning the fastq file for one condition (6h INFa treatment). Generally each ChIP-seq experiment has one file measuing the signal and one file measuring the input (cross-linked and sonicated but not immuno-precipitated).

| STAT1  | control |
| ------------- | ------------- |
| STAT1_30m_IFNa.fastq  | INP_30m_IFNa.fastq  |
| STAT1_6h_IFNa.fastq  | INP_6h_IFNa.fastq  |

The other files have been processed as downsampled BAM files. They have already been aligned and processed. We are looking at STAT1_6h_IFNa.fastq
```
cd ~/Day3
mkdir -p results 
cd raw_data
head STAT1_6h_IFNa.fastq
```
As you can see, we have the fastq standard of four lines per read. These reads are short, so bowtie1 or bwa might be more suitable, but we use bowtie2 because read lengths are generally > 50. 
```
@SRR502327.1 HWI-EAS233:6:1:3:1856 length=27
TTTATCTTGTTNNACCATCCGTACAAT
+SRR502327.1 HWI-EAS233:6:1:3:1856 length=27
BACA=C<,?B(%%<*CB8BCC/C9,A<
```
### Alignment with bowtie2
| arguments  | definition |
| ------------- | ------------- |
| -q  | fastq input file  |
| -p  2 |  2 threads used to do alignment |
| -k  1 | searching for at most one alignment (not optimal)  |
| --local  | local alignment  |
| --no-unal | supress output for records that didn't align | 

```
bowtie2 -q -p 2 -k 1 --local --no-unal -x genome/hg19 raw_data/STAT1_6h_IFNa.fastq > results/STAT1_6h_IFNa.sam
```


Now we have a SAM file and some QC output.

![alt text](../img/alignment_bowtie2.png)

Normally, we would go right ahead and turn it into a bam file, but we need to normalize it. Imagine the sequencing depth of the other experiments is different. If we are looking at the significance of peaks, we need to downsample reads so they are the same between experiments and same between controls. In this case, the STAT1 experiments were downsampled to 11,000,000 reads and the input were downsampled to 19,000,000.

![alt text](../img/sam_output.png)

```
cd results
less STAT1_6h_IFNa.sam # q to exit less 





```
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
