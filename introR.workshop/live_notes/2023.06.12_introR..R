#Packages
#Comments are like this. Put anything you want!

## Install a package from CRAN
# install.packages("tidyverse")
## Load a package
library(tidyverse)

## Install a package from Bioconductor
# install.packages("BiocManager")
# BiocManager::install("limma")
library(limma)

# Set a seed for reproducibility
set.seed(352)

# Loading simple data
meta_csv <- read.table(file="data/meta.csv", sep=",", header=TRUE)
##Get help
?read.table

# A simpler option with smart defaults
meta_csv <- read.csv(file = "data/meta.csv")

# Load complex data
load("data/meta.RData")

# Data types
## Explore your data
## Dimensions
dim(meta)
ncol(meta)
nrow(meta)

## A single variable/column. Called a vector
meta$ptID

## Data types
class(meta$ptID)
class(meta$age_dys)

class(meta$condition)
## Factors have specific, limited levels
levels(meta$condition)
## Factor formating is not saved in a csv/tsv/table. You have to use RData
class(meta_csv$condition)

## More classes
class(meta$RNAseq)
class(meta$total_seq)

## Get all classes for the data frame variables
str(meta)
glimpse(meta)

## Accessing complex data types
## A list of data frames
load("data/RSTR_data_clean_subset.RData")
class(dat)
class(dat$genes)
## Just keep useing $ to go deeper into the list object
dat$genes$hgnc_symbol

#Some data types use @ instead of $

# Working with data
## Use R for math and stats
var(meta$total_seq)

## Ask which samples (rows) have more than 10 million sequences
meta$total_seq > 10E6
## Ask which samples have methylation data
meta$methylation == TRUE

## Get unique values
unique(meta$sex)
## Note how factors also show you the levels
unique(meta$condition)

## Using incorrect data type. It may run without an error, but the result is not what you want
mean(meta$libID)

## Subsetting data
## Access rows in a data frame
### [ rows , columns ]
meta[7, ]
meta[c(7,8), ]

## Filter just Media samples
meta$condition
meta$condition == "Media"
logical.vector <- meta$condition == "Media"
meta[logical.vector, ]

## Filter just Mtb samples
logical.vector2 <- meta$condition == "Mtb"
logical.vector2
meta[logical.vector2, ]

## Access columns of a data frame
meta[ , 3]
meta[, "condition"]

## Working with missing data
# missing data NA is different from the work "NA"

# Plot
## Histogram of total sequences
hist(meta$total_seq)
## Histogram of just Media samples
plot1 <- hist(meta[ logical.vector, ]$total_seq)
plot(plot1)

# The tidyverse way
## Checkout https://r4ds.had.co.nz/ or
## https://bigslu.github.io/workshops/
meta %>% 
  filter(condition == "Media") %>% 
  pull(total_seq) %>% 
  hist()

# Saving data
## One data frame
meta_media <- meta[ logical.vector, ]
write.csv(x = meta_media, file = "data/meta_media.csv")

## Many data objects
save(meta, meta_media, plot1, 
     file = "data/meta_all.RData")

## Save a PNG plot
## measures in pixels
png(filename="plot1.png", width = 500, height=250) 
hist(meta$total_seq)
dev.off()
