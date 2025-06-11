#### Install packages ####
#CRAN packages
install.packages(c("tidyverse", "ggrepel", "scales"), 
                 Ncpus=4)

#Bioconductor packages
install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "biomaRt", "patchwork"))

#GitHub packages
install.packages("devtools")
devtools::install_github("zhangyuqing/sva-devel")
devtools::install_github("BIGslu/kimma")
devtools::install_github("BIGslu/BIGpicture")
devtools::install_github("BIGslu/RNAetc")

#### Load packages ####
library(tidyverse)
library(scales)

#### Load data ####
flagstat <- kimma::example.seasnake$flagstat
picard <- kimma::example.seasnake$picard
# patient <- kimma::example.seasnake$patient

patient <- kimma::example.seasnake$patient %>% 
  mutate(asthma = fct_relevel(factor(asthma), "healthy", after = 0)) 
sample <- kimma::example.seasnake$sample %>% 
  mutate(virus = fct_relevel(factor(virus), "none", after = 0)) %>% 
  mutate(batch = factor(batch))

#merge metadata
meta <- full_join(sample, patient, by = "ptID") %>% 
  full_join(flagstat, by = "libID") %>% 
  full_join(picard, by = "libID")

#counts data
count <- kimma::example.seasnake$fcounts

#Non-kimma data
# count <- read_tsv("path/user name/count/fcounts.tsv")

#### Quality control ####
#Set cutoffs
seq_cutoff <- 1E6
cv_cutoff <- 1
align_cutoff <- 75

#total sequences plot
ggplot(meta, aes(x = reorder(libID, QC_pass), 
                 y = QC_pass)) +
  geom_col() #geometry


ggplot(meta, aes(x = reorder(libID, QC_pass), y = QC_pass)) +
  geom_col() +
  #Add cutoff line
  geom_hline(yintercept = seq_cutoff) +
  #Log scale y-axis
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  #Beautify
  theme_classic() +
  labs(x = "Library", y = "Pass-filter sequences (log scale)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
