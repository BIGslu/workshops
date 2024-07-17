# Comments are like this! Not code. Put anything here

# Install packages from CRAN
install.packages("tidyverse", Ncpus = 4)
install.packages("BiocManager")

# Install packages from Bioconductor
BiocManager::install("limma")

# Load packages
## Has to be done EVERY time you open RStudio
library(tidyverse)
library(limma)

# Load data
## Tabular data (csv, tsv, etc)
meta_csv <- read.table(file = "data/meta.csv", sep=",", header=TRUE)
### A better way since we know it's a csv
meta_csv <- read.csv(file = "data/meta.csv")
## Tab delimited sep="\t"

## RData
### Positives: can be multiple objects, is compressed, keeps the object names, can be non-tabular data, and more!
load("data/meta.RData")

# Help function
?read.table

# Explore data
dim(meta)
class(meta)

## Extract a column
meta$ptID
## Another way
# meta[rows, columns]
meta[ , "ptID"]

## Types of data
class(meta$ptID) #character
class(meta$age_dys) #integer
class(meta$total_seq) #numeric
class(meta$condition) #factor
levels(meta$condition) ##factor order
class(meta$RNAseq) #logical i.e. TRUE/FALSE

## Look at all the data
str(meta)

# Working with data
## Statistics
var(meta$total_seq)
mean(meta$age_dys)
mean(meta$libID) #warning: incorrect class

## Questions about the data
### Which rows have greater than 10 million sequences
meta$total_seq > 10E6
### What are the unique values
unique(meta$sex)
unique(meta$condition)

## Changing a data class
as.factor(meta$sex)
### Save change to data frame. Lots of as.XXX() functions
meta$sex <- as.factor(meta$sex)

## Subsetting data
### Row 2, column 3. Best to worst
meta$condition[2]
meta[2, "condition"]
meta[2, 3]

### only media samples
meta[c(1,3,5,7,9,11,13,15,17,19), ] #terrible!

logical_vector <- meta$condition == "Media"
logical_vector
meta[logical_vector, ] #much better!

### only media samples with > 10E6 sequences
logical_vector2 <- meta$condition == "Media" & meta$total_seq > 10E6
logical_vector2
meta[logical_vector2, ]

# Plotting
## base R
hist(meta$total_seq)
boxplot(meta$condition, meta$total_seq)

## ggplot
ggplot(data = meta, aes(x=condition, y=total_seq)) +
  geom_boxplot()

### Add themes for a quick beautification!
ggplot(data = meta, aes(x=condition, y=total_seq)) +
  geom_boxplot() +
  theme_bw()

ggplot(data = meta, aes(x=condition, y=total_seq)) +
  geom_boxplot() +
  theme_classic()

ggplot(data = meta, aes(x=condition, y=total_seq)) +
  geom_boxplot() +
  theme_dark()
