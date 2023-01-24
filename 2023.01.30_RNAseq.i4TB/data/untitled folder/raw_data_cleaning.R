#### Setup ####
library(tidyverse)
library(stringr)

meta <- read_csv("0_data/2020.11.20RSTR_Hawn_metadata.csv")
count <- read_csv("0_data/Hawn_UgandaMacrophage_counts.csv") %>% 
  rename("hgnc_symbol"=`...1`)

#### Format data ####
#List samples with RNAseq data
samp <- data.frame(libID_orig=colnames(count)) %>% 
  separate(libID_orig, into=c("RS_SUB_ACCESSION_NO","condition"), sep="_",
           remove=FALSE) %>% 
  mutate(condition = recode(condition, "TB"="Mtb","MEDIA"="Media"))

meta.ltbi <- meta %>% 
  #Filter LTBI samples with RNAseq data
  filter(KCHCA_AGE_YR_CURRENT >= 18 & Sample_Group == "LTBI" & 
           RS_SUB_ACCESSION_NO %in% samp$RS_SUB_ACCESSION_NO) %>% 
  #Select variables of interest
  select(FULLIDNO, RS_SUB_ACCESSION_NO, M0_KCVSEX,
         KCHCA_AGE_YR_CURRENT, mono.RNAseq, methylation) %>%
  #Convert T/F variables
  mutate(mono.RNAseq = ifelse(is.na(mono.RNAseq),FALSE,TRUE),
         methylation = ifelse(is.na(methylation),FALSE, TRUE)) %>% 
  #Add library ID for MEDIA and TB
  inner_join(samp) %>% 
  #rename to short names
  rename(sex=M0_KCVSEX, age_yrs=KCHCA_AGE_YR_CURRENT,
         RNAseq=mono.RNAseq) 
  
# Select 10 random donors
set.seed(8759)
samp10 <- sample(unique(meta.ltbi$RS_SUB_ACCESSION_NO), size=10, 
                 replace=FALSE)

meta.ltbi10 <- meta.ltbi %>% 
  filter(RS_SUB_ACCESSION_NO %in% samp10) %>% 
  #De-identify
  mutate(ptID = paste0("pt", 
                       str_pad(as.numeric(factor(RS_SUB_ACCESSION_NO)),2,pad = "0")),
         libID = paste(ptID, condition, sep="_")) %>% 
  #add variables for introR
  mutate(age_dys = age_yrs*365,
         ptID_old = paste0("pt", 
                           str_pad(as.numeric(factor(RS_SUB_ACCESSION_NO)),5,pad = "0")))
  

#Filter and rename counts table
count.ltbi <- count %>% 
  pivot_longer(-hgnc_symbol, names_to = "libID_orig") %>% 
  inner_join(select(meta.ltbi10, libID_orig, libID)) %>% 
  select(-libID_orig) %>% 
  pivot_wider(names_from = "libID")

#Add total sequences
tot.seqs <- count.ltbi %>% 
  pivot_longer(-hgnc_symbol, names_to = "libID") %>% 
  group_by(libID) %>% 
  summarise(total_seq = sum(value, na.rm=TRUE))

#### Kinship ####
kin <- read_csv("0_data/kinship_Hawn_all.csv") %>% 
  inner_join(distinct(meta.ltbi10, FULLIDNO, ptID), 
             by=c("rowname"="FULLIDNO")) %>% 
  select(-rowname) %>% 
  rename(rowname=ptID) %>% 
  pivot_longer(-rowname) %>% 
  inner_join(distinct(meta.ltbi10, FULLIDNO, ptID), 
             by=c("name"="FULLIDNO")) %>% 
  select(-name) %>% 
  arrange(ptID) %>% 
  pivot_wider(names_from = "ptID") %>% 
  arrange(rowname)

#### Save ####
write_csv(count.ltbi, file="0_data/raw_counts.csv")

meta.ltbi10 %>% 
  left_join(tot.seqs) %>% 
  select(libID, ptID, condition, age_dys, sex, ptID_old, RNAseq,
         methylation, total_seq) %>% 
  arrange(ptID, condition) %>% 
  write_csv(file="0_data/metadata.csv")

write_csv(kin, file="0_data/kinship.csv")

