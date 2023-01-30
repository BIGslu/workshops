#Packages
library(tidyverse)
library(limma)
library(kimma)

#### Data ####
dat <- example.voom

#Explore data
dat$E[1:3,1:3]
dat$targets[1:3,]
dat$genes[1:3,1:4]
dat$weights[1:3,1:3] #from voomWithQualityWeights

#### Limma: simple linear model ####
#Build model
mm_limma <- model.matrix(~virus + asthma, data = dat$targets)
head(mm_limma)

#Fit model
fit_limma <- lmFit(object = dat$E, design = mm_limma,
                   weights = dat$weights)

#Estimate p-values
efit_limma <- eBayes(fit = fit_limma)

#View top results
fdr_limma <- topTable(fit = efit_limma)
View(fdr_limma)

#View results for individual variables (ie coefficients)
##Virus
fdr_limma_hrv <- topTable(fit = efit_limma, number = Inf,
                      coef = "virusHRV")
View(fdr_limma_hrv)

#How many genes significant for virus at FDR < 0.05
fdr_limma_hrv %>% 
  filter(adj.P.Val < 0.05) %>% 
  nrow()

##Asthma
fdr_limma_asthma <- topTable(fit = efit_limma, number = Inf,
                          coef = "asthmaasthma")
View(fdr_limma_asthma)

##All variables in one data frame
fdr_limma <- extract_lmFit(design = mm_limma, fit = efit_limma)

####limma: paired design ####
#Calculate gene expression correlation between paired samples
consensus.corr <- duplicateCorrelation(object = dat$E,
                                       design = mm_limma,
                                       block = dat$targets$ptID)
consensus.corr$consensus.correlation

#Use paired correlation in model fit
fit_limma2 <- lmFit(object = dat$E, design = mm_limma,
                   weights = dat$weights,
                   block = dat$targets$ptID,
          correlation = consensus.corr$consensus.correlation)
#Estimate p-values
efit_limma2 <- eBayes(fit = fit_limma2)

#View results
fdr_limma2 <- extract_lmFit(design = mm_limma, 
                            fit = efit_limma2)
View(fdr_limma2)

####limma: compare models ####
#Format results for merging
temp <- fdr_limma$lm %>% 
  select(gene, variable, FDR) %>% 
  filter(variable == "virusHRV") %>% 
  rename(`limma` = FDR)

temp2 <- fdr_limma2$lm %>% 
  select(gene, variable, FDR) %>% 
  filter(variable == "virusHRV") %>% 
  rename(limma_dupCor = FDR)

#Merge results and plot FDR
full_join(temp, temp2) %>% 
  ggplot(aes(x = `limma`, y = limma_dupCor)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_classic() +
  coord_fixed() +
  labs(title = "FDR for HRV vs media")

#plot sigma
##lower sigma is better model fit
data.frame(`limma` = fdr_limma$lm.fit$sigma,
           limma_dupCor = fdr_limma2$lm.fit$sigma) %>% 
  
  ggplot(aes(x = `limma`, y = limma_dupCor)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_classic() +
  coord_fixed() +
  labs(title = "Sigma")

#### kimma: simple and linear mixed effects models ####
fit_kimma <- kmFit(dat = dat, 
                   model = "~ virus + asthma + (1|ptID)", 
                   use.weights = TRUE,
                   run.lm = TRUE, run.lme = TRUE,
                   metrics = TRUE)
## No need for Bayes or extract_fit. Results are immediately available
View(fit_kimma)
fit_kimma$lm.fit %>% head()
fit_kimma$lm %>% head()

#### kimma: compare models ####
fit_kimma_all <- full_join(fit_kimma$lm.fit, 
                           fit_kimma$lme.fit, 
                           by = c("gene"), 
                           suffix = c("_lm","_lme"))

fit_kimma_all %>%
  ggplot(aes(x = AIC_lm, y = AIC_lme)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_classic() +
  coord_fixed() +
  labs(title = "AIC")

head(fit_kimma_all)

#How many gene better fit by simple lm
fit_kimma_all %>% 
  filter(AIC_lm < AIC_lme) %>% 
  nrow()

#Mean AIC across genes
mean(fit_kimma$lm.fit$AIC)
mean(fit_kimma$lme.fit$AIC)

#Cumulative AIC of all genes
sum(fit_kimma$lm.fit$AIC)
sum(fit_kimma$lme.fit$AIC)

#How many genes are significant for each variable?
fit_kimma$lme %>% 
  filter(FDR < 0.05) %>% 
  count(variable)

summarise_kmFit(fit_kimma$lme)

#### Pathways ####
library(SEARchways)
library(BIGpicture)

#### DEG enrichment ####
#List signif genes
genes.virus <- fit_kimma$lme %>% 
  filter(variable == "virus" & FDR < 0.1) %>% 
  pull(gene)
genes.virus

#run enrichment on Hallmark database
enrich.virus <- BIGprofiler(
  gene_list = list("virus"=genes.virus),
  ID = "ENSEMBL", category = "H")
head(enrich.virus)

#plot signif enrichment
plot_enrich(enrich.virus, fdr.cutoff = 0.5,
            fdr.colors = c(0.05,0.5))

#run enrichment on GO terms
enrich.virus2 <- BIGprofiler(
  gene_list = list("virus"=genes.virus),
  ID = "ENSEMBL", category = "C5", subcategory = "GO:BP")

#plot signif enrichment
plot_enrich(enrich.virus2, fdr.cutoff = 0.2,
            fdr.colors = c(0.05,0.2))

#### GSEA ####
#Extract estimates (ie fold change) for each gene for virus variable
FC.virus <- fit_kimma$lme %>% 
  filter(variable == "virus") %>% 
  select(variable, gene, estimate)
View(FC.virus)

#Set seed to reproducibly deal with ties
set.seed(3524786)
#Run GSEA on Hallmark database
gsea.virus <- BIGsea(gene_df = FC.virus,
                     ID = "ENSEMBL",
                     category = "H")

head(gsea.virus)

#Plot significant pathways
p1 <- plot_gsea(gsea.virus, fdr.cutoff = 0.2, 
          fdr.colors = c(0.05,0.2))
class(p1)
p1

#You can continue to edit the plots if desired
p1 + theme_dark()

#### Custom gene sets ####
#db	= If not using Broad databases, a data frame with gene ontology including gene set name (column 1: gs_name) and gene ID (column2: gene_symbol, entrez_gene, or ensembl_gene as matches your gene_list names)

#Format a custom database
set.seed(78)
db_custom <- data.frame(
  gs_name = c('gs1','gs2'),
  ENSEMBL = c(
    paste(sample(dat$genes$geneName, 10), collapse=","),
    paste(sample(dat$genes$geneName, 20), collapse=",")
  )
) %>% 
  #Format to lists
  rowwise() %>% 
  mutate(ENSEMBL = str_split(ENSEMBL, ","))

#Another way to format a custom database
set.seed(78)
db_custom2 <- data.frame(
  gs_name = c('gs1','gs2'),
  ENSEMBL = sample(dat$genes$geneName, 50)
  ) %>% 
  group_by(gs_name) %>% 
  summarise(ENSEMBL = list(ENSEMBL))

head(db_custom)

#Run GSEA
gsea.custom <- BIGsea(gene_df = FC.virus,
                     ID = "ENSEMBL",
                     db = db_custom)
#Plot results
plot_gsea(gsea.custom, fdr.cutoff = 1)

#### Fin ####
