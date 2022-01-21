#2021-10-28 Intro R

# Install packages from CRAN
install.packages(pkgs="ggplot2")
?install.packages

#First try installing NOT from source
#Then try from source with compilation

#To run code
#Mac - Cmd-Shift-Enter
#PC - Crtl-Shift-Enter

#Load packages
library(package=ggplot2)
library(ggplot2)
?library

# Install packages from Bioconductor
install.packages("BiocManager")
BiocManager::install("limma")

#results in the same thing as above
library(BiocManager)
install("limma")

#Load package. Should see no messages if successful
library(limma)

#In general
### Update all (a) if asked
### But still don't compile from source

#Example update message when installing
# Seeing this: The downloaded binary packages are in
# /var/folders/yp/vtvlf1gs23l89sqy8p5blp800000gp/T//RtmppPirgx/downloaded_packages
# Old packages: 'lattice', 'mgcv', 'nlme', 'survival'
# Update all/some/none? [a/s/n]: 

# Potential errors
# Look for whichever package is the issue (here rlang) and then install that separately before retrying install of the original package

# The downloaded binary packages are in
# C:\Users\user\AppData\Local\Temp\RtmpmgtcZb\downloaded_packages
# Warning message:
#   In file.copy(savedcopy, lib, recursive = TRUE) :
#   problem copying C:\Users\user\Documents\R\win-library\4.1\00LOCK\rlang\libs\x64\rlang.dll to C:\Users\user\Documents\R\win-library\4.1\rlang\libs\x64\rlang.dll: Permission denied
install.packages("rlang")

#Load data
#csv with column names in first row
meta <- read.table(file="data/RSTR_meta_subset.csv", header=TRUE, sep=",")
#See help page for how to use tsv, etc
?read.table

#Load RData
load(file="data/RSTR_data_clean_subset.RData")

#Working with data
#What type of data?
class(dat)
class(meta)
#dimensions?
dim(meta)
dim(dat)

#get a column from data frame
meta$libID
#Get a column from a complex list
dat$targets$libID

#Look at classes again
#surprise! targets in dat is the same data frame as the meta we loaded from csv
class(meta$libID)
class(meta$lib.size)
class(dat$targets$libID)

#other classes
##factor
##integer
##logical

#Calculate basic stats
median(meta$norm.factors)
#unique values/levels
unique(meta$condition)

#class type has to be correct for a function to work
#example mean of a character does not make sense
mean(meta$libID)

#Ask R questions
#where is lib.size greater than 1E7
meta$lib.size > 1E7
#where is it less than 1E7 AND greater than 1E8 (answer: none because that's impossible)
meta$lib.size < 1E7 & meta$lib.size > 1E8


#basic format: meta[rows, columns]
# is blank, no filtering done
#another way to get a column
meta[,"lib.size"]

