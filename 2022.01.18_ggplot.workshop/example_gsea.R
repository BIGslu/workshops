library(tidyverse)
#Broad gene set database
library(msigdbr)
#GSEA
library(fgsea)

#Data
library(kimma)
dat <- example.voom

# Broad database
H <- as.data.frame(msigdbr(species = "Homo sapiens", 
                           category = "H")) %>% 
  #Keep gene ID that match expression data gene ID
  select(gs_name, ensembl_gene) %>% 
  #Collapse all genes in each gene set into 1 row each
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(ensembl_gene))) %>%
  #Convert to list
  deframe()

# Calculate fold change
#Extract expression data
FC <- as.data.frame(dat$E) %>% 
  #Move gene IDs from rownames to a column
  rownames_to_column("ensembl_gene_id") %>% 
  #Make long format with all expression values in 1 column
  pivot_longer(-ensembl_gene_id, 
               names_to = "libID", values_to = "expression") %>% 
  #Add metadata
  left_join(dat$targets) %>% 
  select(ensembl_gene_id, expression, donorID, virus) %>% 
  #Make wide with media and tb expression in separate columns
  pivot_wider(names_from = virus, values_from = expression) %>% 
  #Calculate tb minus media fold change (delta for change)
  #Because expression values are log2, subtraction is the same as division
  mutate(delta = HRV-none) %>% 
  #Calculate mean fold change per gene
  group_by(ensembl_gene_id) %>% 
  summarise(mean.delta = mean(delta, na.rm=TRUE)) %>% 
  #Arrange by descending fold change
  arrange(desc(mean.delta))

#Vector of mean fold change values
FC.vec <- FC$mean.delta
#Add names
names(FC.vec) <- FC$ensembl_gene_id

# Run GSEA
gsea.H <- as.data.frame(fgseaSimple(pathways = H, 
                                    stats = FC.vec,
                                    nperm = 1000,
                                    scoreType = "std"))

write_csv(gsea.H, "data/example_gsea.csv")
