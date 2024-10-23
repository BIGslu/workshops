library(tidyverse)

# Load data
# Similar to base R but with some bells and whistles
# meta <- read.csv(...)
meta <- read_csv(file="data/metadata.csv")
class(meta)

# Package: dplyr
# Select a column in base R
meta$libID
meta[,"libID"]
meta[,1] #Danger!

# Select two columns in dplyr
select(meta, libID, condition)

# Remove a column in dplyr
select(meta, -libID)

# Select everything else
## Useful in reordering columns
## This is select ptID and then everything else in the original order
select(meta, ptID, dplyr::everything())
## After everything(), you can also remove a column
select(meta, ptID, dplyr::everything(), -ptID_old)

# Filter rows
## Media samples
filter(meta, condition=="Media")
## Media samples with > 10 million sequences
filter(meta, condition=="Media" & total_seq>10E6)

# rename columns
## format is always new name = old name
rename(meta, tb=condition)

# String things together, the bad way2
## This creates a cluttered environment
meta_clean <- select(meta, libID, condition)
meta_clean2 <- filter(meta_clean, condition=="Media")
## And this makes in unclear what meta_clean contains 
meta_clean <- select(meta, libID, condition)
meta_clean <- filter(meta_clean, condition=="Media")

# Pipes %>% allow you to string tidyverse functions together without saving the intermediate data frames
meta_clean <- select(meta, libID, condition) %>% 
  filter(condition=="Media")
# If works by taking the output of the first function (a data frame) and using that
# as the input (data) of the next function
# This is why all tidyverse functions have the first parameter as data

## select(meta, libID, condition) %>% 
##   filter("assumed data is here", condition=="Media")

# Like to always use pipes, even for the first function
meta_clean <- meta %>% 
  # Select columns of interest
  select(libID, condition) %>% 
  # Filter media samples
  filter(condition=="Media")

# base R pipe /> is also can option and something you may see out in the wild

# renaming bad column names examples
## By hand
rename(meta, pTNF_cd14=`%TNF+CD14+`)
## Or try janitor's default renaming
library(janitor)
meta_rename <- clean_names(meta)
colnames(meta_rename)

# Calculate summary statistics
## Can do with base R
mean(meta[meta$condition=="Media",]$total_seq)
mean(meta[meta$condition=="Mtb",]$total_seq)

# With dplyr, can do for each group all at once
meta %>% 
  group_by(condition) %>% 
  summarise(mean_seq = mean(total_seq, na.rm=TRUE),
            n = n())

# Dealing with NA types
# Remember to look at your data and make sure NA are coded correctly
meta <- read_csv(file="data/metadata.csv",
                 na=c(" ",NaN, NA))
#Looking at example NA column
# R recognized NA will be grey italics in the preview
meta %>% 
  mutate(test=NA) %>% View

#Infinity in R
Inf
-Inf
## Suggestion filter to remove Inf or -Inf before summarise otherwise the results
## won't be meaningful

# More on groups
# Ungroup the data. Until you ungroup, the data are assumed to remain grouped
# which can impact downstream analyses
## Option 1
meta %>% 
  group_by(condition) %>% 
  summarise(mean_seq = mean(total_seq)) %>% 
  ungroup()
## Option 2 currently being assessed by Posit for it's usefulness
meta %>% 
  group_by(condition) %>% 
  summarise(mean_seq = mean(total_seq),
            .groups = "drop")

# create new column aka mutate
mutate(meta, age_yrs = age_dys/365) %>% 
  View()

# Package: tidyr
# Make a wide table (fewer rows, more columns)
meta_wide <- meta %>% 
  select(-libID) %>% 
  pivot_wider(names_from=condition, 
              values_from=total_seq)
View(meta_wide)

# Make a long table (more rows, fewer columns)
meta_long <- meta_wide %>% 
  pivot_longer(c(Media, Mtb),
               names_to = "condition",
               values_to = "total_seq")
View(meta_long)

# Naming pivot so the new columns are more clear
meta %>% 
  select(-libID) %>% 
  pivot_wider(names_from=condition, 
              values_from=total_seq,
              names_prefix="total_seq_") %>% 
  View()

# More dplyr
# Joining tables
# https://stat545.com/join-cheatsheet.html
# BiocManager::install("limma")

# Load data with multiple data frames
load("data/dat_voom.RData")
str(dat)
# metadata
View(dat$targets)
# gene expression data
View(dat$E)

as.data.frame(dat$E) %>% 
  #move rownames into a data column so they are not lost
  rownames_to_column("hgnc_symbol") %>% 
  #pivot everything except hgnc_symbol
  pivot_longer(-hgnc_symbol) %>% 
  #Join sample metadata
  full_join(dat$targets, by=c("name"="libID")) %>% 
  #Select 1 gene of interest
  filter(hgnc_symbol=="TNF") %>% 
  
  #plot
  ggplot(aes(x=condition, y=value)) +
  geom_boxplot()

# Dealing with dates
library(lubridate)
## Note you do not have these data so I am providing the outputs in comments
class(meta$VISIT)
# [1] "character"
meta$VISIT[1]
# [1] "9/19/17"

meta <- meta %>% 
  #Format the date
  mutate(visit_date = as.Date(VISIT, "%m/%d/%y")) %>% 
  select(subjid, VISIT, visit_date)
class(meta$visit_date)
# [1] "Date"
meta$visit_date[1]
# [1] "2017-09-19"

# Use lubridate to easily do math with dates
# Days between dates
days(meta$visit_date[1]-meta$visit_date[2])
# [1] "172d 0H 0M 0S"

# Days in the month of this date
days_in_month(meta$visit_date[1])
# Sep 
# 30

# grepl example
# Select all unformatted dates that start with 7/
# this is a regex (regular expression)
meta %>% filter(grepl("^7\\/", VISIT))
