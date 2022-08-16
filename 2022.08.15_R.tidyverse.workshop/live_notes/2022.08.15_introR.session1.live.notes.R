# Install packages
## Only need to do once on your computer. Commented out (#) here so does not run

# install.packages("tidyverse")
# install.packages("ggrepel")
# install.packages("BiocManager")
# BiocManager::install("limma")

# Load packages
## Need to complete every time you open R/RStudio
library(tidyverse)
library(ggrepel)
library(limma)

# Getting help
## For a specific function name
?average
?mean
## Search functions and help pages for key word(s)
??average

# Load data
meta <- read.table(file = "data/metadata.csv", 
                   header = TRUE, sep = ",")
## Note could also use read.csv which assumes correct header and sep parameters

# Explore the data
## Table dimensions
dim(meta)
## Pull out a single column by its name
meta$ptID

## Types of data in R
class(meta$ptID)
class(meta)
class(meta$age_dys)
class(meta$total_seq)
class(meta$RNAseq)
meta$RNAseq

#See all info on a data frame
str(meta)

# Complex data type
## Load
load("data/dat_voom.RData")
class(dat)

## Pull out nested pieces with $
dat$targets
dat$targets$libID

## These data are S3. You may encounter S4 in R, which uses @ instead of $

# Some simple math and stats
mean(meta$total_seq)
sd(meta$total_seq)
meta$total_seq > 10E6
table(meta$total_seq > 10E6)

unique(meta$condition)
unique(dat$targets$condition)

# Subsetting data
## [rows, columns]

meta[1,1]
meta[meta$total_seq > 10E6, ]

## Save the TRUE/FALSE vector to use in subsetting. 
## Outcome is the same as the code above but it's easiter to read when you have
## complete statements
logical.vector <- meta$total_seq > 10E6
logical.vector
meta[logical.vector, ]
