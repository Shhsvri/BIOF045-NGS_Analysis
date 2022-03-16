# This is intro to R
# 07/12/2021
# Shahin Shahsavari

# I am adding 3 and 5
3 + 5

# These are numeric values
a <- 3
a = 3

a + 5

# These are characters
c = "CCR6"

character_a = "3"

# These are logical

logical_a = TRUE

# Numeric Vectors

vector_a = c(3, 4, 6, 10)

# Character Vectors

vector_b = c("a", "b","c", "d")
vector_c = c("a", "b", "c", 1)


# Data Frame

df_a = data.frame(vector_a, vector_b, vector_c)
View(df_a)

## Extract columns from a dataframe

df_a$vector_a
mode(df_a$vector_a)
mode(df_a$vector_b)

## subset a vector

vector_a[3]
vector_b[2]

## subset a dataframe
df_a[2]
df_a[1, 3]
View(df_a)

## Find dimentions of a dataframe
ncol(df_a)
nrow(df_a)
dim(df_a)

## Find dimentions of a vector
length(vector_a)
length(df_a$vector_a)

## column names and row names
colnames(df_a)
rownames(df_a)

## How to assign names to rows in a dataframe

names = c("Sample1", "Sample2", "Sample3", "Sample4")
rownames(df_a) = names
View(df_a)


# using test dataframe called mtcars

View(mtcars)
dim(mtcars)
nrow(mtcars)
ncol(mtcars)
mtcars[2,4]
mtcars$hp[2]

## Access values using column names and row names
mtcars["Mazda RX4 Wag", "hp"]


# Functions in R

numeric_vector = c(45, 53, 22, 110, 55)
sum(numeric_vector)
mean(numeric_vector)
sqrt(81)
sqrt(numeric_vector)

sum(sqrt(numeric_vector))

## write your own function in R

square_it <- function(x) {
  square = x * x
  return(square)
}

square_it(21)

# Packages in R
## sometimes we need extra functions developed by the community
## We install them this way

install.packages("PACKAGE_NAME")
install.packages("readxl")

library(tidyverse)

# conditional filtering
filter(mtcars, gear == 4)

# subsetting a dataframe
select(mtcars, hp)

# Pipes %>%
sqrt(numeric_vector) %>% sum()

# head function with pipe
filter(mtcars, gear == 4) %>% head()
filter(mtcars, gear == 4) %>% head(n = 9)

# filter() then select()

filter(mtcars, hp > 140) %>% select(mpg, cyl, wt, gear, carb)

# extract all but one column
filter(mtcars, hp > 140) %>% select(-vs)

# Reading data into R from a file

setwd("~/Day1/R_data")
getwd()
list.files()

# read.table( file = "PATH")
counts = read.table(file = "counts.csv", sep = ",")
View(counts)

dim(counts)
colnames(counts)

# look at the first column, and take the sum of all reads
sum(counts$sample2)

# take the mean of RPKMs in the first column: "sample2"
mean(counts$sample2)

# filter out the counts that are less than 10, and equal to 0 for sample2

filter(counts, sample2 > 10) %>% nrow()

## != means does not equal
filter(counts, sample2 != 0) %>% nrow()



# Read the metadata file
metadata = read.table(file = "mouse_exp_design.csv", sep = ",")
View(metadata)

## extract the rows where genotype is equal to "Wt"
filter(metadata, genotype == "Wt")


## extract the rows where genotype is "Wt" and celltype is typeA
genotype == "Wt") %>% filter(celltype == "typeA")
typeA_Wt = filter(metadata, genotype == "Wt" & celltype == "typeA")

View(typeA_Wt)

# write.table
write.table(typeA_Wt, file = "metadata_typeA_Wt.csv", sep = ",")
