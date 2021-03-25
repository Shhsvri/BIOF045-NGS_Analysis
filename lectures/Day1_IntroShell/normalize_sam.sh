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