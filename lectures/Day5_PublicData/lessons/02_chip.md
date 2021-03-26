# ChIP-seq analysis 
We are working with the same data as before. 

### Alignment (DO NOT RUN)
We are performing a local alignment to get better contiguous alignments between our reads and the genome

Arguments for bowtie2
```
-q (fastq input)
-p 2 (multiprocessing)
-k 1 find first alignment (suboptimal but fast)
--no-unal ignore unaligned reads
-x genome/hg19 genome location
-X The maximum fragment length for valid paired-end alignments (not included)
```
Code 
```Bash
bowtie2 -q -p 2 -k 1 --local --no-unal -x genome/hg19 raw_data/STAT1_6h_IFNa.fastq > results/STAT1_6h_IFNa.sam
```

### Read normalization (DO NOT RUN)
When comparing peak sets, we want the comparison to be fair, to have started from a shared baseline. This is why it's important to downsample reads. 
The reads in the counts for the TFs must be the same as each other, and the reads in the input controls must be the same.
![alt text](../img/mappable.png)

Duplicates only need to be removed if there is a lot of low quality reads (nonredundant fraction is high).  Otherwise, you may be removing some of your signal.

I showed you a script you could write for normalizing reads. The script takes in two arguments, the SAM file and the number of reads to downsample to. 
```Bash
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
After saving as `normalize_sam.sh`, the script could be run like the following:
```Bash
./normalize_sam.sh STAT1_30m_IFNa.sam 11000000
```
Note that permissions sometimes have to be changed to run scripts; `chmod 777 normalize_sam` should work. The above command would downsample STAT1_30m_IFNa.sam to 11M reads. 

These downsampled SAM files could then be sorted and indexed in another script:
```Bash
# get the base name of the SAM file
baseSam=basename $1 .sam

# convert to BAM, sort, and index
samtools view -b $1 > tmp
samtools sort tmp > $baseSam.bam
samtools index $baseSam.bam
```
After saving this as `samToBam.sh`, it could be run like this:
```
./samToBam norm_STAT1_30m_IFNa.sam
```
Remember that permissions need to be changed in some cases. 

### Peak calling (Run again)
Finally we covered peak calling with MACS3, which uses multiple Possion models to model counts across the genome.

Arguments for MACS3
```
# -t treatment
# -c control
# -n name prefix
# --outdir output directory
# -g human sample
# --bdg generate bedGraph
```

Commands for MACS
```Bash
macs3 callpeak -t norm_STAT1_30m_IFNa.bam -c norm_INP_30m_IFNa.bam -n STAT1_30m_IFNa --outdir ../peaks -g hs --bdg -q 0.05 -f BAM
macs3 callpeak -t norm_STAT1_6h_IFNa.bam -c norm_INP_6h_IFNa.bam -n STAT1_6h_IFNa --outdir ../peaks -g hs --bdg -q 0.05 -f BAM
```

