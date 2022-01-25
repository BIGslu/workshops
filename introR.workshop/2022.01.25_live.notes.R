#############
# Introduction to R
# 2022.01.25
# https://github.com/BIGslu/workshops/tree/main/introR.workshop
# https://youtu.be/pCotIcLfet4
############
# Style guide
# https://style.tidyverse.org/

# Comments start with a #
# line 2 if you need one

# Basic R syntax
# my_function = function(data = ..., output = ...)

# R packages (CRAN)
## Install. Once per computer
install.packages("ggplot2")
## Load package
library(ggplot2)

# R package (Bioconductor)
install.packages("BiocManager") # Package for install Bioconductor packages
BiocManager::install("limma")
##Or
library(BiocManager)
install("limma")

# R is case sensitive
## These are two different objects in the Environment
a <- 1
A <- 2
## <- is the preferred assignment operator but = also works
a = 1
A = 2

# Help function
## by function name
?install.packages
?average #error
## more general search when you don't know the function name
??average

# Load data
## No formatting
read.table(file = "data/RSTR_meta_subset.csv")
## Correct csv with column names formatting
read.table(file = "data/RSTR_meta_subset.csv", sep = ",", header = TRUE)

## Parameter order is important if you don't use the parameter names.
## If order if different from defaults (see help page), the function will not work
read.table(file = "data/RSTR_meta_subset.csv", ",", TRUE) # error
## But if order is the same as defaults, it works
## Caution when doing something like this as others (and probably future you)
## won't know exactly what this does without consulting the help page.
read.table("data/RSTR_meta_subset.csv", TRUE, ",")
## One except is that the first parameter is usually data or a file name
## so removing its parameter name is common
read.table('data/RSTR_meta_subset.csv', sep = ',', header = TRUE)

## Save the table to an object in the environment so we can use it!
meta <- read.table("data/RSTR_meta_subset.csv", sep = ",", header = TRUE)

# Another type of data
## R formatted data
## RData can contain non-data frames, multiple objects, and are compressed
## to reduce space usage
load("data/meta.RData")
load("data/RSTR_data_clean_subset.RData")

# Data types: standard
class(meta)
## Access a column in a data frame
class(meta$lib.size)
class(meta$libID)

## more case sensitive notes
## To be a logical class, TRUE and FALSE must be used
True #not correct
TRUE

# Data types: complex 
## S3
## Access with $
class(dat)
class(dat$targets)
class(dat$targets$lib.size)

## S4
## Use @ instead of $

# Do math in R
## See stats package for more basic functions
mean(meta$lib.size)
var(meta$lib.size)
## Ask T/F question
meta$lib.size > 10E6
## Other useful questions
unique(meta$condition)
## or if it were a factor
levels(meta$condition)

# Subsetting
## The hard way based on position. Not reproducible and very error prone
meta$lib.size[c(4,5,10,11,13,14,15,17,18,19,20)]
meta$lib.size[4]

## The better way using a T/F vector
logical.vector <- meta$lib.size > 10E6
logical.vector
## access rows with the format meta[rows, columns]
meta[logical.vector, ]

# Other operators (also see notes)
# < > <= >=
# == equal
# != not equal
# & and
# | or

## String together multiple T/F statements
logical.vector <- meta$lib.size > 10E6 & meta$condition == "MEDIA"
meta[logical.vector, ]

## Use the subset data in another function
meta.media <- meta[meta$condition == "MEDIA", ]
mean(meta.media$lib.size)
