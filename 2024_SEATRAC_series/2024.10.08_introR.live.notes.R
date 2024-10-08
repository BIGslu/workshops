#### Intro to R ####
# SEATRAC workshop series
# 2024-10-08

#### Load packages ####
# Unless you tell R which package to use, it can't find packages that are not in base R
separate()
# Note that tidyverse is the a package
tidyverse::separate()
# tidyr (which is part of tidyverse) is the package for this function
tidyr::separate()
# But this is annoying so you can load all functions in a package to your workspace with library()
library(tidyverse)
separate()
# Another package loading
library(limma)

#### Load data - csv ####

read.table(file = "data/meta.csv", sep = ",",
           header = TRUE)
#Not in an Rproj? Try this file path for file =
# "~/Desktop/introR/data/meta.csv"

# Parameter order matters if you don't provide parameter names like file = .
## This is a wrong order
read.table("data/meta.csv", ",", TRUE)
# This is the right (default) order
read.table("data/meta.csv", TRUE, ",")
# How I actually do it. Often people skip the name for the first parameter as it is usually the data (file or object in your environment)
read.table("data/meta.csv", sep = ",",
           header = TRUE)
# But this just prints to the console below

#Save to environment (so you can do things with the data!)
meta_csv <- read.table("data/meta.csv", 
                       sep = ",",
                       header = TRUE)

#### Getting help ####
# When you know the function name
?read.table
?separate
# When you don't know the function name and just want to serve all the help pages
??read.table

#### Load data - RData ####
load("data/meta.RData")
# Similar for Rds

#### Removing data ####
# Sweeper button in the Environment tab 
# Or using code
rm(meta_csv)

#### Data types - simple ####
# Info on a data frame
dim(meta)
class(meta)

# Get a column from a data frame
meta$ptID
class(meta$ptID)

# Factors
class(meta$condition)
levels(meta$condition)
meta$condition

# Numbers
class(meta$age_dys)
class(meta$total_seq)

# Logical (true or false)
class(meta$RNAseq)
meta$RNAseq
## R accepts
TRUE
T
FALSE
F
# R doesn't accept
True
true

# Info on an entire data frame
str(meta)

#### Data types - complex ####
# S3 list
load("data/RSTR_data_clean_subset.RData")
class(dat)
# Get 1 data frame from the list
class(dat$genes)
head(dat$genes)

# Get 1 column from 1 data frame from the list
dat$genes$geneName

# All data frames in the list
names(dat)

#### Working with data ####
# Run basic statistics
mean(meta$total_seq)
var(meta$total_seq)
sd(meta$total_seq)
median(meta$total_seq)
# etc...

#### Subsetting rows ####
# Get media samples
## [rows, columns]
## BAD - using indexing
meta_media <- meta[c(1,3,5,7,9,11,13,15,17,19), ]
## BETTER - using a conditional statement
meta$condition == "Media"
logical_vector <- meta$condition == "Media"
meta[logical_vector, ]

# media samples with > 10 million seqs
logical_vector2 <- meta$condition == "Media" & meta$total_seq > 1E7
meta[logical_vector2, ]

#### Subsetting columns ####
# Similar to rows but in second slot in [ , ]
# meta[, put_something_here]
meta[, 2]

#### Plotting in base R ####
# Let R guess the plot type
plot(meta$total_seq)
# Histogram
hist(meta$total_seq)

#### Plotting in ggplot ####
# Teaser of ggplot! Stay tuned for more workshops

ggplot(meta, aes(x=condition, 
                 y=total_seq, 
                 color=condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() 

# Different themes (i.e. different looks)
ggplot(meta, aes(x=condition, 
                 y=total_seq, 
                 color=condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw()

ggplot(meta, aes(x=condition, 
                 y=total_seq, 
                 color=condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_minimal()
