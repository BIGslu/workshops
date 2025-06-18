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
library(ggrepel)
# library(biomaRt)
library(patchwork)
library(edgeR)
library(limma)
#
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

#Median CV coverage
#variability of gene coverage

#Set cutoffs
CV_max <- max(1, max(meta$MEDIAN_CV_COVERAGE))
cv_cutoff_label <- 0.9
align_cutoff_label <- 75

ggplot(meta, aes(x = MEDIAN_CV_COVERAGE, y = paired_mapped/QC_pass*100)) +
  geom_point() +
  #Rescale axis limits
  lims(x = c(0,CV_max), y=c(0,100)) +
  #Add cutoff lines
  geom_hline(yintercept = align_cutoff, lty="dashed") +
  geom_vline(xintercept = cv_cutoff, lty="dashed") +
  #Label points outside cutoffs
  geom_text_repel(data=filter(meta, MEDIAN_CV_COVERAGE > cv_cutoff_label | 
                                paired_mapped/QC_pass*100 < align_cutoff_label),
                  aes(label=libID), show.legend = FALSE, max.overlaps = Inf,
                  box.padding = 1) +
  #Beautify
  theme_classic() +
  labs(x = "Median CV coverage", y="Percent alignment")

#remove poor quality libraries
meta.filter <- meta %>% 
  filter(MEDIAN_CV_COVERAGE < cv_cutoff & QC_pass > seq_cutoff &
           paired_mapped/QC_pass*100 > align_cutoff)
#nothing was removed b/c all high-quality
nrow(meta)
nrow(meta.filter)

count.filter <- count %>% 
  select(1, all_of(meta.filter$libID))

#### Protein-coding ####
#Get database
ensembl <- biomaRt::useEnsembl(biomart="ensembl", 
                               dataset="hsapiens_gene_ensembl")

#Format gene key
key <- biomaRt::getBM(attributes=c("ensembl_gene_id", 
                                   "hgnc_symbol",
                                   "gene_biotype", "chromosome_name",
                                   "start_position", "end_position"),
                      mart=ensembl) %>% 
  #Filter protein coding genes
  filter(gene_biotype == "protein_coding")

key.filter <- key %>% 
  #Filter protein coding genes in count table
  filter(ensembl_gene_id %in% count$Geneid) %>% 
  #collapse multiannotations.
  group_by(ensembl_gene_id, gene_biotype, 
           chromosome_name, start_position, end_position) %>% 
  summarise(symbol = list(unique(hgnc_symbol)), .groups = "drop")

#Filter pc genes
count.filter.pc <- count.filter %>% 
  filter(Geneid %in% key.filter$ensembl_gene_id)

dim(count.filter.pc)
dim(meta.filter)

#### Batch effects ####
BIGpicture::plot_pca(count.filter.pc, meta=meta.filter, 
                     vars="batch", transform_logCPM=TRUE)
#with scaling. recommended 
BIGpicture::plot_pca(count.filter.pc, meta=meta.filter, 
                     vars="batch", transform_logCPM=TRUE,
                     scale=TRUE)

#label duplicates
unique.ID <- c("ptID", "virus")

#Identify duplicate libraries
dups <- meta.filter %>% 
  unite("dupID", all_of(unique.ID), remove=FALSE)  %>% 
  count(dupID) %>% 
  filter(n>1)

#Create duplicate ID variable
meta.filter.dups <- meta.filter %>% 
  unite("dupID", all_of(unique.ID),remove=FALSE)  %>% 
  mutate(duplicate = ifelse(dupID %in% dups$dupID, dupID, NA))

BIGpicture::plot_pca(count.filter.pc, meta=meta.filter.dups, 
                     vars="duplicate", transform_logCPM=TRUE,
                     scale=TRUE)
# Combat-seq
count.filter.pc.combat <- count.filter.pc %>% 
  #Convert count df to matrix
  column_to_rownames("Geneid") %>% 
  as.matrix() %>% 
  #Batch correction
  sva::ComBat_seq(., batch = meta.filter$batch,
                  group = meta.filter$virus,
                  covar_mod = model.matrix(~ sex, meta.filter),
                  #don't use on real data. just to make run fast
                  shrink = TRUE, gene.subset.n = 2) %>% 
  as.data.frame() %>% 
  rownames_to_column("Geneid")

BIGpicture::plot_pca(count.filter.pc.combat, 
                     meta=meta.filter.dups, 
                     vars=c("batch","duplicate"),
                     transform_logCPM=TRUE, scale=TRUE) %>% 
  wrap_plots(ncol=2)

#### PCA outliers ####
pca3a <- BIGpicture::plot_pca(count.filter.pc.combat, 
                              meta=meta.filter, 
                              vars="outlier", transform_logCPM=TRUE,
                              scale = TRUE)
pca3a
BIGpicture::plot_pca(count.filter.pc.combat, 
                     meta=meta.filter, 
                     vars="outlier", transform_logCPM=TRUE,
                     scale = TRUE, outlier_sd = 1.5)
#no need to filter any outliers
not.outlier <- pca3a$outlier$data %>% 
  filter(col.group == "no")

meta.filter.out <- meta.filter %>% 
  filter(libID %in% not.outlier$libID)

count.filter.pc.combat.out <- count.filter.pc.combat %>% 
  select(1, all_of(meta.filter.out$libID))

#### Remove duplicates ####
#Find libraries with highest seqs
meta.filter.out.dedup <- meta.filter.out %>% 
  group_by_at(unique.ID) %>% 
  slice_max(order_by = QC_pass)

#Filter lower seq libraries
count.filter.pc.combat.out.dedup <- count.filter.pc.combat.out %>% 
  select(1, all_of(meta.filter.out.dedup$libID))
dim(count.filter.pc.combat.out.dedup)

#### normalization ####
#Order metadata by library ID
meta.filter.out.dedup.ord <- meta.filter.out.dedup %>% 
  arrange(libID)

#Order counts by library ID and gene ID
count.filter.pc.combat.out.dedup.ord <- count.filter.pc.combat.out.dedup %>% 
  select(1, all_of(meta.filter.out.dedup.ord$libID)) %>% 
  arrange(Geneid)

#Order gene key by gene ID
key.filter.ord <- key.filter %>% 
  arrange(ensembl_gene_id)

#check libraries
identical(meta.filter.out.dedup.ord$libID,
          colnames(count.filter.pc.combat.out.dedup.ord)[-1])
#Check genes
identical(key.filter.ord$ensembl_gene_id,
          count.filter.pc.combat.out.dedup.ord$Geneid)

# Move gene names from a variable in the df to rownames and format to matrix
count.filter.pc.combat.out.dedup.ord.mat <- 
  count.filter.pc.combat.out.dedup.ord %>% 
  column_to_rownames("Geneid") %>% 
  as.matrix()

dat <- DGEList(
  #count table
  counts=count.filter.pc.combat.out.dedup.ord.mat,
  #metadata
  samples=meta.filter.out.dedup.ord,
  #gene info
  genes=key.filter.ord)

# filter low abundance genes
BIGpicture::plot_mv(dat, design = "~ virus")
dat.abund <- RNAetc::filter_rare(dat, min_CPM = 0.1, min.sample = 3,
                                 gene.var="ensembl_gene_id")
BIGpicture::plot_mv(dat.abund, design = "~ virus")

#too harsh
dat.abund2 <- RNAetc::filter_rare(dat, min_CPM = 1, min.sample = 3,
                                 gene.var="ensembl_gene_id")
BIGpicture::plot_mv(dat.abund2, design = "~ virus")

#not harsh enough
dat.abund2 <- RNAetc::filter_rare(dat, min_CPM = 0.1, min.sample = 1,
                                  gene.var="ensembl_gene_id")
BIGpicture::plot_mv(dat.abund2, design = "~ virus")

#look at removed genes
rare <- as.data.frame(cpm(dat$counts)) %>% 
  #Filter genes removed
  rownames_to_column("ensembl_gene_id") %>% 
  filter(!(ensembl_gene_id %in% rownames(dat.abund$counts))) %>% 
  #Add gene symbols
  left_join(dat$genes, by = "ensembl_gene_id") %>% 
  #format
  select(-c(chromosome_name:end_position)) %>%
  pivot_longer(-c(ensembl_gene_id, gene_biotype, symbol)) %>% 
  group_by(ensembl_gene_id, gene_biotype, symbol) %>% 
  summarise(mean.CPM = mean(value, na.rm=TRUE),
            min.CPM = min(value, na.rm=TRUE),
            max.CPM = max(value, na.rm=TRUE),
            express.in.libs = length(value[value > 0]),
            .groups="drop")

dim(dat.abund)

write_csv(rare, file="data_clean/rare_genes.csv")
rare

# TMM trimmed mean of means
dat.abund.norm <- calcNormFactors(dat.abund, method = "TMM")
View(dat.abund.norm$counts)
View(dat.abund.norm$samples)
dat.abund.norm$samples$norm.factors

dat.abund.norm.voom <- voom(
  dat.abund.norm,
  design=model.matrix(~ virus, 
                      data=dat.abund.norm$samples),
  plot=TRUE)
View(dat.abund.norm.voom$E)
names(dat.abund.norm.voom)

dat.abund.norm.voom <- voomWithQualityWeights(
  dat.abund.norm,
  design=model.matrix(~ virus, 
                      data=dat.abund.norm$samples),
  plot=TRUE)
View(dat.abund.norm.voom$weights)

#### limma batch effects ####
#See notes https://bigslu.github.io/tutorials/RNAseq/2.RNAseq_counts.to.voom.html#limma_batch2

#### Final PCA check ####
BIGpicture::plot_pca(dat.abund.norm.voom, 
                              vars="outlier",
                              scale = TRUE)
BIGpicture::plot_pca(dat.abund.norm.voom, 
                     vars=c("virus","sex",
                            "ptID","age"),
                     scale = TRUE) %>% 
  wrap_plots(ncol=2)
