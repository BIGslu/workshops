# Install packages
## Only need to do once on your computer. Commented out (#) here so does not run

## CRAN
# install.packages("tidyverse")
# install.packages("ggrepel")

## Bioconductor
# install.packages("BiocManager")
# BiocManager::install("limma")

# Load packages
library(tidyverse)
library(ggrepel)
library(limma)

# Load data frame
meta <- read.table(file = "data/metadata.csv", 
                   header = TRUE, sep = ",")
## You can remove parameter names if you know they are in the correct, default order
meta <- read.table("data/metadata.csv", 
                   TRUE, ",")
## I tend to youse parameter names expect for the first input, which is usually
## data or a file
meta <- read.table("data/metadata.csv", 
                   header = TRUE, sep = ",")

## Open the table
View(meta)

# Help
?read.table

# Load RData
load("data/dat_voom.RData")

# Explore data
## Dimensions
## rows columns
dim(meta)
dim(dat)

## From data frames
class(meta)
## You extract columns with $
meta$ptID

# From S3 data type
class(dat)
## You also use $ but it can have 2+ levels of nesting
dat$targets$ptID

# Data types
class(meta$ptID)
class(meta$age_dys)
class(meta$total_seq)
class(meta$RNAseq)

## See all info for data in a table
str(meta)

# Working with and subsetting
## Some math and stats
meta$total_seq
mean(meta$total_seq)
sd(meta$total_seq)

## Ask R for a TRUE/FALSE vector for a statement
meta$total_seq > 10E6

## Use that T/F vector to filter rows that are TRUE
logical.vector <- meta$total_seq > 10E6
#[ rows , columns ]
meta[logical.vector, ]

## Further filter to just get 1 column by name
meta[logical.vector, "sex"]
