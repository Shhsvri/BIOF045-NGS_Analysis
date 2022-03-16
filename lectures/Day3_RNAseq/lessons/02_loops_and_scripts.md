---
title: "The Shell: Loops & Scripts"
author: "Shahin Shahsavari"
date: "November 2021"
---

Approximate time: 60 minutes

## Learning Objectives

* Capture previous commands into a script to re-run in one single command
* Understanding variables and storing information
* Learn how to use variables to operate on multiple files


Now that you've been introduced to a number of commands to interrogate your data, wouldn't it be great if you could do this for each set of data that comes in, without having to manually re-type the commands?

Welcome to the beauty and purpose of shell scripts.

## Shell scripts

Shell scripts are **text files that contain commands we want to run**. As with any file, you can give a shell script any name and usually have the extension `.sh`. For historical reasons, a bunch of commands saved in a file is usually called a shell script, but make no mistake, this is actually a small program. 


We are finally ready to see what makes the shell such a powerful programming environment. We are going to take the commands we repeat frequently and save them into a file so that we can **re-run all those operations** again later by typing **one single command**. Let's write a shell script that will do two things:

1. Tell us our current working directory
2. List the contents of the directory 

Before we do so, let's try a new command for *printing* text into the terminal:

```bash
$ echo "Hello World"
```

> `echo` displays a line of text onto your terminal. You can redirect these into files.

First open a new file using `vim`:

```bash
$ cd ~/Day3/scripting
$ vim listing.sh
```

Then type in the following lines in the `listing.sh` file:

```bash
#!/bin/bash

# This is my first shell script
# Shahin Shahsavari
# 03/16/2022

echo "Your current working directory is:"
pwd
echo "These are the contents of this directory:"
ls -l 
```

> **NOTE** lines that start with <kbd>#</kbd> are usually comments. The one exception to this rule is the shebang *#!* which tells our
system how to interpret a script. The other lines are ignored and treated as commends.
It is crucial to document what your scripts are supposed to do so you and your
colleagues could read them. This helps make your code reproducible.

Exit `vim` and save the file. Now let's run the new script we have created. To run a shell script you usually use the `bash` or `sh` command.

```bash
$ source listing.sh
```

> Did it work like you expected?
> 
> Were the `echo` commands helpful in letting you know what came next?

This is a very simple shell script, just to introduce you to the concept. Later in this session, we will be learning how to write more complex scripts to illustrate the power of scripting and how it can make our lives (when coding) much easier. Any type of data you will want to analyze will inevitably involve not just one step, but many steps and perhaps many different tools/software programs. Compiling these into a shell script is the first step in creating your analysis workflow!

Before we jump into more scripts, we will take a moment to cover some key concepts to help get you there.

## Bash variables
A *variable* is a common concept shared by many programming languages. Variables are essentially a symbolic/temporary name for, or a reference to, some information. Variables are analogous to "buckets", where information can be stored, maintained and modified without too much hassle. 

Extending the bucket analogy: the bucket has a name associated with it, i.e. the name of the variable, and when referring to the information in the bucket, we use the name of the bucket, and do not directly refer to the actual data stored in it.

Let's start with a simple variable that has a single number stored in it:

```bash
$ num=25
```

*How do we know that we actually created the bash variable?* We can use the `echo` command to print to terminal:

```bash
$ echo num
```

What do you see in the terminal? The `echo` utility takes what arguments you provide and prints to terminal. In this case it interpreted `num` as a a character string and simply printed it back to us. This is because **when trying to retrieve the value stored in the variable, we explicitly use a `$` in front of it**:

```bash
$ echo $num
```

Now you should see the number 25 returned to you. Did you notice that when we created the variable we just typed in the variable name? This is standard shell notation (syntax) for defining and using variables. When defining the variable (i.e. setting the value) you can just type it as is, but when **retrieving the value of a variable don't forget the `$`!** 

Variables can also store a string of character values. In the example below, we define a variable or a 'bucket' called `file`. We will put a filename `ptA.fastq` as the value inside the bucket.

```bash
$ file=ptA.fastq
```

Once you press return, you should be back at the command prompt. Let's check what's stored inside `file`, but first move into the `raw_fastq` directory::

```bash
$ echo $file
```

Let's try another command using the variable that we have created. We can also count the number of lines in `Mov10_oe_1.subset.fq` by referencing the `file` variable:

```bash
$ wc -l $file
```

> *NOTE:* The variables we create in a session are system-wide, and independent of where you are in the filesystem. This is why we can reference it from any directory. However, it is only available for your current session. If you exit the cluster and login again at a later time, the variables you have created will no longer exist.

---

**Exercise**

1. change the file variable and set it equal to `ptB.fastq`
2. reuse the file variable and print the number of lines in the fastq file
---

## Loops

Another powerful concept in the Unix shell and useful when writing scripts is the concept of "Loops". We have just shown you that you can run a single command on multiple files by creating a variable whose values are the filenames that you wish to work on. But what if you want to **run a sequence of multiple commands, on multiple files**? This is where loop come in handy!

Looping is a concept shared by several programming languages, and its implementation in bash is very similar to other languages. 

The structure or the syntax of (*for*) loops in bash is as follows:

```bash
for (variable_name) in (list)
do
	(command $variable_name)
	.
	.
done
```

where the ***variable_name*** defines (or initializes) a variable that takes the value of every member of the specified ***list*** one at a time. At each iteration, the loop retrieves the value stored in the variable (which is a member of the input list) and runs through the commands indicated between the `do` and `done` one at a time. *This syntax/structure is virtually set in stone.* 


#### What does this loop do? 

```bash
cd ~/Day3/scripting
for x in *.fastq
 do
   echo $x
   ls -l $x
 done
```

We will use for loops and shell scripts to align our RNA seq samples tomorrow.

For the DNA data, we could use the following shell script to iterate over both samples and generate the vcf files.

```bash
#!/bin/bash

# BIOF045: 03/16/2022
# This script for DNA alignment, sorting, and indexing

## 0. set up the file structure change your directory

cd ~/Day3/scripting
genome=~/Day2/genome/hg38.fa


for i in *.fastq
 do
	file=${i%.*}
	## 1. BWA MEM alignment with 2 threads
	
	bwa mem -t 2 \
		$genome \
		$i > $file.sam
	
	
	## 2. Convert sam to bam
	
	samtools view -b $file.sam > $file.bam
	
	
	## 3. Sort your bam file using samtools
	
	samtools sort $file.bam > $file.sorted.bam
	
	
	## 4. markduplicates with PICARD
	
	cd results
	
	PicardCommandLine MarkDuplicates \
		I=$file.sorted.bam \
		O=$file.markdup.sorted.bam \
		M=$file_md_metrics.txt
	
	
	## 5. Index the bam file
	###	after this step you could view the bam file in IGV
	
	samtools index $file.markdup.sorted.bam
	
	
	## 6. Generate the VCF file using bcftools
	
	bcftools mpileup -f $genome $file.markdup.sorted.bam | bcftools call -mv -Ov -o $file.vcf
	
	## 7. Remove the unneeded files that were generated during the alignment
	
	rm $file.sam $file.bam $file.sorted.bam
 done

```

---
