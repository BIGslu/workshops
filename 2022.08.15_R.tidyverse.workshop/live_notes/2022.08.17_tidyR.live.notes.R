# Load packages
library(tidyverse)

# Load data
meta <- read_csv(file = "data/metadata.csv")
class(meta)
View(meta)

# dplyr
## select columns
select(meta, ptID) # A single column
select(meta, ptID, total_seq) # 2 columns
select(meta, -ptID_old) # Remove 1 column
select(meta, -c(ptID_old, sex)) # Remove multiple columns
select(meta, -ptID_old, -sex) #Equivalent way to remove multiple columns

## filter rows
filter(meta, condition == "Media") # Keep media libaries only
meta[meta$condition == "Media", ] # Base R way to do the same thing

filter(meta, condition == "Media" & condition == "Mtb") #no rows returned b/c condition cannot be BOTH (&) Media and Mtb for a single library
filter(meta, condition == "Media" | condition == "Mtb") #all rows returned b/c all libraries are either (|) Media or Mtb 

## From a question
## How to filter for all libraries where "Media" exists in the condition value
## Would keep conditions like Media1, Media2, etc
filter(meta, grepl("Media", condition))
## Check out regular expressions (regex) for ways to search for basically can pattern of letters, numbers, and/or symbols

## rename a column
### new name = old name
rename(meta, infection = condition)

## summarize / summarise 
## all rows mean of age
summarise(meta, mean(age_dys))
## all rows mean and std deviation of age
summarise(meta, 
          mean_age_dys = mean(age_dys),
          sd_age_dys = sd(age_dys))

## group_by
## Does not change the meta data frame is how it looks
meta_grouped <- group_by(meta, condition)
## But under the hood, R knows the groups
class(meta_grouped)
## summarise mean age within each condition group
summarise(meta_grouped, mean_age_dys = mean(age_dys))

## pipes %>% (tidyverse) or /> (base R)
#cmd+shift+m shortcut to type a pipe
### For those familiar with math notation
### f( ) %>% g( ) = f(g( ))

### or more specially, the pipe places the resulting data frame as the first argument in the next function
### Here the . is where the pipe places it
group_by(meta, condition) %>% 
  summarise(., mean_age_dys = mean(age_dys))

### This is the same without a pipe
summarise(group_by(meta, condition),
          mean_age_dys = mean(age_dys))

## mutate new variables
### age in years. New column appears on the far right of the data frame
mutate(meta, age_yrs = age_dys/365.25)

### here is the same thing again with a pipe
meta <- meta %>% 
  mutate(age_yrs = age_dys/365.25)

## Exercises
#1. Create a new column for total sequences in million total_seq_millions based on total_seq.
meta %>% mutate(total_seq_million = total_seq/1000000)
# Or use scientific notation for large numbers
meta %>% mutate(total_seq_million = total_seq/1E6)

#2. What is another way to end up with only the Media rows instead of condition == "Media"?
filter(meta, condition != "Media")
filter(meta, condition %in% c("Mtb"))

#3. Try calculating the mean total number of sequences for the Media and Mtb conditions.
summarise(meta_grouped, 
          mean_total_seq = mean(total_seq))

## Or from meta before we grouped it in meta_grouped
meta %>% 
  group_by(condition) %>% 
  summarise(mean_total_seq = mean(total_seq)) %>% 
  ungroup()

# tidyr
## Widen data
### Take the unique levels in the name variable (condition) and give them each a new column
### Fill those columns with the data in the value variable (total_seq)
meta_wide <- pivot_wider(meta, 
                         names_from = condition,
                         values_from = total_seq)
## Oh no! We still have two rows per patient because libID is unique for each libary

## Remove unique variables so can put all libraries from 1 patient in 1 row
meta_wide <- meta %>% 
  select(-libID) %>% 
  pivot_wider(names_from = condition,
              values_from = total_seq)

## Optionally add to the column names to better understand what the data are
meta_wide <- meta %>% 
  select(-libID) %>% 
  pivot_wider(names_from = condition,
              values_from = total_seq,
              names_prefix = "total_seq_")

## See we now have half the rows
dim(meta)
dim(meta_wide)

## Undo the pivot back to long format
## Notice we have to tell it which columns to merge back into just 1 name and 1 value column
meta_long <- meta_wide %>% 
  pivot_longer(c(total_seq_Media, total_seq_Mtb),
               names_to = "condition",
               values_to = "total_seq")
dim(meta_long) # 1 column missing b/c we removed libID before
dim(meta)

## Bringing it all together
## Let's show a bunch of what we've learned in 1 piped string
meta_clean <- meta %>% 
  select(-libID) %>% 
  filter(sex == "F") %>% 
  rename(infection = condition) %>% 
  mutate(age_yrs = age_dys/365.25) %>% 
  pivot_wider(names_from = infection,
              values_from = total_seq)

# Back to dplyr
## Load expressiondata
load("data/dat_voom.RData")
## This is a list of data frames
names(dat)

dat$E[1:3,1:3] # log2 counts per million
dat$targets # Roughly the same as the meta data frame we've been using thus far

## Note that the expression data are in a matrix with row names
class(dat$E)
rownames(dat$E)

## Tidyverse needs a data frame without row names
E_long <- as.data.frame(dat$E) %>% 
  rownames_to_column("gene") %>% 
  ## Get libID in a column for matching to the metadata
  pivot_longer(-gene, names_to = "libID",
               values_to = "log2CPM")

## Note how this changed the data
dim(dat$E)
dim(E_long)

## Now both the data frames have a column for libID
E_long[1:3,]
dat$targets[1:3,]

## merge data frames
## if you don't give column(s) to merge by, joins use all columns with matching names
full_data <- inner_join(dat$targets, E_long) 
## Or you can specify the column(s) to match rows by
full_data <- inner_join(dat$targets, E_long,
                        by = "libID")
## If your columns did not have the same name, you can specify
## Will not run. Just an example of syntax
inner_join(dat$targets, E_long,
           by = c("left_name"="right_name"))

# Exercises
#1. Filter the metadata to just Media samples. Without running it, outline what the data would look like if you inner_join the long expression data E_long with metadata containing only Media samples.
media_data <- dat$targets %>% 
  filter(condition == "Media") %>% 
  inner_join(E_long, by="libID")
dim(full_data)
## Notice fewer rows b/c Mtb library expression is removed
dim(media_data)
table(is.na(media_data))

#2. Similarly, what would happen with full_join?
media_data2 <- dat$targets %>% 
  filter(condition == "Media") %>% 
  full_join(E_long, by="libID")
dim(full_data)
dim(media_data2)
## notice rows with NA where metadata are missing
table(is.na(media_data2))
