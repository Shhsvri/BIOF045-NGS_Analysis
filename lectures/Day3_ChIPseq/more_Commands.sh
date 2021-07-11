# grab the data 
mkdir raw_data
cd raw_data

fasterq-dump -Z SRR502329 -e 16
fasterq-dump -Z SRR502327 -e 16
fasterq-dump -Z SRR502228 -e 16
fasterq-dump -Z SRR502225 -e 16

mv SRR502329.fastq STAT1_30m_IFNa.fastq
mv SRR502327.fastq STAT1_6h_IFNa.fastq
mv SRR502228.fastq INP_30m_IFNa.fastq
mv SRR502225.fastq INP_6h_IFNa.fastq

cd ..
mkdir sam_data 

wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
unzip hg19.zip

# -q fastq 
# -p 4 threads 
# -k search for at most this many alignments 
# --local local alignment 
# --no-unal ignore unaligned SAM records
# -x basename of genome 
bowtie2 -q -p 4 -k 1 --local --no-unal -x hg19 raw_data/STAT1_30m_IFNa.fastq > sam_data/STAT1_30m_IFNa.sam
bowtie2 -q -p 4 -k 1 --local --no-unal -x hg19 raw_data/STAT1_6h_IFNa.fastq > sam_data/STAT1_6h_IFNa.sam
bowtie2 -q -p 4 -k 1 --local --no-unal -x hg19 raw_data/INP_30m_IFNa.fastq > sam_data/INP_30m_IFNa.sam
bowtie2 -q -p 4 -k 1 --local --no-unal -x hg19 raw_data/INP_6h_IFNa.fastq > sam_data/INP_6h_IFNa.sam

./normalize_sam.sh STAT1_30m_IFNa.sam 11000000
./normalize_sam.sh STAT1_6h_IFNa.sam 11000000
./normalize_sam.sh INP_30m_IFNa.sam 19000000
./normalize_sam.sh INP_6h_IFNa.sam 19000000

curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3 get-pip.py
pip install numpy --upgrade --ignore-installed
pip install macs3 --upgrade --ignore-installed

bash samToBam.sh norm_STAT1_30m_IFNa.sam
bash samToBam.sh norm_STAT1_6h_IFNa.sam 
bash samToBam.sh norm_INP_30m_IFNa.sam
bash samToBam.sh norm_INP_6h_IFNa.sam

# -t treatment
# -c control
# -n name prefix
# --outdir output directory
# -g human sample 
macs3 callpeak -t norm_STAT1_30m_IFNa.bam -c norm_INP_30m_IFNa.bam -n STAT1_30m_IFNa --outdir ../peaks -g hs --bdg -q 0.05 -f BAM
macs3 callpeak -t norm_STAT1_6h_IFNa.bam -c norm_INP_6h_IFNa.bam -n STAT1_6h_IFNa --outdir ../peaks -g hs --bdg -q 0.05 -f BAM