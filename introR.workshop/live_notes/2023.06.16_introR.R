# Packages
# Notes to self, others
## Install packages from CRAN
## Only need to run once of your computer
install.packages("tidyverse") 
## Need to run every time you open R
library(tidyverse)

# Install from Bioconductor
install.packages("BiocManager")
BiocManager::install("limma")
# Answer "a" to updating and "no" to compilation
library(limma)

# Set seed for reproducibility
set.seed(345869)

# Load tabular data
meta_csv <- read.table(file = "data/meta.csv", sep = ",", header = TRUE)

## Open a help page
?read.table

## Use logical defaults
meta_csv <- read.csv(file = "data/meta.csv")

# Read in RData
load(file = "data/meta.RData")

# Explore the data
dim(meta) # Lists rows and then columns
nrow(meta)
ncol(meta)

# Access a variable
meta$ptID

# Types of data
class(meta)
class(meta$ptID)
class(meta$condition)
## more about factors
levels(meta$condition)
meta$condition
## csv does not have factor formatting
class(meta_csv$condition)

## more data types
class(meta$age_dys)
class(meta$total_seq)
class(meta$RNAseq) # TRUE vs FALSE

# see all the data
str(meta)
glimpse(meta)

# Working with data
## Running basic stats
mean(meta$total_seq)
sd(meta$total_seq)
var(meta$total_seq)
## More help pages when you don't know the function name
??average

## Ask questions
## Rows with more than 10 million sequences
meta$total_seq > 10E6
## Unique values of characters or factors
unique(meta$sex)
unique(meta$condition)

# Using indexes of data
meta$libID[1]
meta[1,1] #rows, columns
meta[1,"libID"]

#large data, see just the start
head(meta)
meta[1:3,1:3]

# Subset data
logical.vector <- meta$total_seq > 10E6
logical.vector

meta[logical.vector, ] #rows, columns
meta_deep <- meta[logical.vector, ]

# Plot
hist(meta$total_seq)
?hist
hist(meta$total_seq, breaks=10)

ggplot(data = meta) +
  geom_histogram() +
  aes(total_seq) +
  theme_classic() +
  labs(title ="Histogram of total sequences")

# Saving
write.csv(meta_deep, file="data/meta_deep.csv")

png(filename = "total_seq.png", width = 500, height=500)
  hist(meta$total_seq, breaks=10)
dev.off()

## In ggplot (tidyverse), use ggsave()
