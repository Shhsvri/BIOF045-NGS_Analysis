# Software Installation

### How to install bwa and STAR or other bioinformatics softwares
```
cd ~/Day5/installation
```

Most of the bioinformatics softwares are open source, and their source code could be accessed on Github. For example, if we search for `STAR`, we find this [page](https://github.com/alexdobin/STAR).

There, you find information on how to install STAR from scratch

```
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz
# Decompress the file using tar
tar -xzf 2.7.8a.tar.gz
# Move into the folder you just downloaded
cd STAR-2.7.8a
```

Depending on the programming language used to create the software, you need a different approach. But with STAR, bwa and most other programs written in `C` and `C++`, you need to run the following:

```
make STAR
```

If you would like to run this from the current folder, you would use `./STAR` from the current directory.

By default, bash will only find the softwares available in the $PATH variable however.

```
echo $PATH
```

If you want to have STAR be availabe from any directory, you need to add its installation folder to $PATH as follows:

```
PATH=~/Day5/installation/STAR-2.7.8a:$PATH
```


### Installing bwa

```
cd ~/Day5/installation
git clone https://github.com/lh3/bwa.git
cd bwa
make
PATH=~/Day5/installation:$PATH
```



### Installing HOMER

Lastly, we would like to install a software we need for ChIPseq called HOMER

```
mkdir ~/homer
cd ~/homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install hg38
```

Lastly, add this to your $PATH variable

```
PATH=~/homer/bin:$PATH
```














