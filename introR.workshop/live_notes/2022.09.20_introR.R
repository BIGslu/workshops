# Introduction to R
# R-Ladies workshop, 2022.09.20
# Full notes https://bigslu.github.io/workshops/introR.workshop/introR.html

## R packages
### Intall from CRAN
install.packages("tidyverse")
### Install from Bioconductor
install.packages("BiocManager")
BiocManager::install("limma")

## Load packages
### You must run *every* time you open R/RStudio
library(tidyverse)
library(limma)

## Load data
meta_csv <- read.table(file = "data/meta.csv", sep=",", header=TRUE)
meta_csv <- read.csv(file = "data/meta.csv")

## Help
?read.table

## Data types
### See notes for more examples
class(meta_csv)
dim(meta_csv)

meta_csv$ptID
### Character: non-numeric letters and characters
class(meta_csv$ptID)
### Factor: character variable with discrete levels
factor(meta_csv$condition)

### RData
load("data/meta.RData")
### Notice how the RData version saved the factor formatting while csv did not
meta_csv$condition
meta$condition

## Doing math
var(meta$total_seq)
mean(meta$total_seq)
mean(meta$sex) #gives warning

meta$total_seq > 10E6

## Subsetting data
logical_vector <- meta$total_seq > 10E6
logical_vector

### [ rows , columns]
### If one is blank, returns all present in the original data
meta[logical_vector, ]
meta[1,1]
