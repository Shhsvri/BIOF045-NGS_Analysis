---
title: R and Tidyverse
author: Shahin Shahsavari
date: October 2020
duration: 90 minutes
---

# R and Tidyverse

## Functions in R

A function is a set of statements organized together to perform a specific task. R has a 
large number of in-built (nase) functions and the user can create their own functions. You
can also access the functions other people have written from CRAN and Bioconductor.

### built\_in functions in R

```R
# create a sequence of numbers from 32 to 44
seq(32,44)

# add numbers from 41 to 67
sum(41:69)

# find the mean of 22 to 115
mean(22:115)

# create a _vector_ that contains numbers 1, 2, 5, 14
c(1,2,5,14)
```

### Writing your own functions

```R
# Create a function to print squares of numbers in sequence.
new.function <- function(a) {
   for(i in 1:a) {
      b <- i^2
      print(b)
   }
}
```

```R

# Call the function new.function supplying 6 as an argument.
new.function(6)
```

### library functions

In bash, we used a variety of programs such as ls, bwa, cat, samtools, and bcftools.
Those are open source libraries that were developed in the C programming language.
We did not have to write any of the code ourselves.

There is a similar concept in R. If we want to analyze our RNAseq data, we can easily download
a set of functions from open source libraries that are uploaded to CRAN and Bioconductor. We will
learn how to install libararies/packages in R in the next section.

## R Libraries

There is a total of 16445 packages/libraries on CRAN.

### installation

In order to install new packages from CRAN, you need to run:

`install.packages("package name")`

A popular package for reading excel files is:

```R
install.packages("readxl")
```

## Tidyverse

Tidyverse is a very popular package in R. It is developed by the same people who maintain RStudio.

In order to load a library and use its functions, you need to:

```R
library("tidyverse")
```

Every time you open R, you need to run this and load the library functions.

### Reading a text file into a data frame using tidyverse

Earlier, we learned how data columns are separated in a sample text file. If you know you
have a comma delimited file, you could easily use `read_csv("/path/to/file")`:

```R
cars <- read_csv("/data/review/cars.csv")
```

For tab separated files, you would use:

RNA\_counts <- read\_tsv("/data/RNAcounts/counts.rpkm.tsv")

### filter()

If you would like to subset certain rows that satisfy a specific condition, you could use the
filter(data, _condition_) function:

```R
filter(RNA_counts,sample12 > 5)
```
> **Note** This function is available through tidyverse. You must use library(tidyverse) before 
running it


### pipes in tidyverse

If you remember, we used `|` to send the output of one command to another command. This saved us a
lot of time and space. The tidyverse package has defiled an operator in R that gets the same job
done `%>%`:

```R
filter(RNA_counts, sample12 > 5) %>%
	filter(sample11 < 6) %>%
	filter(sample8 != 0)
```

**Exercise**

- read the vcf file we generated yesterday into a data frame in R
- find all the rows that have a quality higher than 20

The file is saved in:

**/data/review/08008.dbSNP.annotate.vcf**

