#### Install packages ####
# install.packages("tidyverse")
# BiocManager::install("limma")
# devtools::install_github("BIGslu/kimma")
# May need to set GitHub credentials if you're logged in
# credentials::set_github_pat(force_new = TRUE)

#### Load packages ####
library(tidyverse)
library(limma)
library(kimma)

#### Load data ####
load("data/dat_voom.RData")
class(dat)
names(dat)

#### limma: RNAseq linear models ####
# Define a model
head(dat$targets)
mm_limma <- model.matrix(~ condition + sex,
                         data = dat$targets)
mm_limma

# Redefine a reference level if desired
dat$targets <- dat$targets %>% 
  mutate(sex = factor(sex, levels=c("M","F")))

mm_limma <- model.matrix(~ condition + sex,
                         data = dat$targets)
mm_limma

# Fit model
fit_limma <- lmFit(object = dat$E, 
                   design = mm_limma,
                   weights = dat$weights)
class(fit_limma)

# Estimate p-values
efit_limma <- eBayes(fit = fit_limma)
class(efit_limma)

# Extract results
fdr_limma <- topTable(fit = efit_limma)
View(fdr_limma)

# Extract results for 1 variable
fdr_limma_mtb <- topTable(fit = efit_limma,
                          coef = "conditionMtb",
                          number = Inf)
View(fdr_limma_mtb)

fdr_limma_sex <- topTable(fit = efit_limma,
                          coef = "sexF",
                          number = Inf)
View(fdr_limma_sex)

# Extract all variables
## This is a kimma wrapper function not available in limma
fdr_limma <- extract_lmFit(design = mm_limma,
                           fit = efit_limma)
names(fdr_limma)
## Model fit metrics
View(fdr_limma$lm.fit)
## Model estimates and significance
View(fdr_limma$lm)

#### Paired design in limma ####
# Calculate mean correlation of genes in media vs genes in mtb
consensus.corr <- duplicateCorrelation(
  object = dat$E,
  design = mm_limma,
  block = dat$targets$ptID)
consensus.corr$consensus.correlation

# Add "paired" design to limma by including mean correlation
fit_limma2 <- lmFit(
  object = dat$E,
  design = mm_limma,
  weights = dat$weights,
  block = dat$targets$ptID,
  correlation = consensus.corr$consensus.correlation)

# Calculate p-values
efit_limma2 <- eBayes(fit_limma2)
#Extract results
fdr_limma2 <- extract_lmFit(design = mm_limma,
                            fit = efit_limma2)

# Compare models with and without "paired" corrrelation
## Format each model's results
temp <- fdr_limma$lm %>% 
  select(gene, variable, FDR) %>% 
  filter(variable == "conditionMtb") %>% 
  rename(limma_FDR = FDR)
temp2 <- fdr_limma2$lm %>% 
  select(gene, variable, FDR) %>% 
  filter(variable == "conditionMtb") %>% 
  rename(limma_FDR2 = FDR)
## Combine and plot
### NOTE: In the workshop, I incorrectly interpreted this like it was sigma (model fit). This is, instead, FDR values. See the workshop notes for interpretation
## https://bigslu.github.io/workshops/2024_SEATRAC_series/3_linear_model_rnaseq.html
full_join(temp, temp2) %>% 
  ggplot(aes(x=limma_FDR, y=limma_FDR2)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  coord_fixed()

#### kimma ####
# DO NOT run unless your have a couple minutes to wait
fit_kimma <- kmFit(
  #Give data
  dat = dat,
  #Define the model
  model = "~ condition + sex + (1|ptID)",
  use_weights = TRUE,
  #Run both the simple and paired models
  run_lm = TRUE, run_lme = TRUE,
  #Calculate fit metrics
  metrics = TRUE
)

# load results if did not run above
load("data/kimma_result.RData")
##Format is similar to extract_lmFit
names(fit_kimma)
View(fit_kimma$lm)
View(fit_kimma$lm.fit)

# Compare kimma models
fit_kimma_all <- full_join(
  fit_kimma$lm.fit,
  fit_kimma$lme.fit,
  by = "gene",
  #Nice trick to have join rename columns for you
  suffix = c("_lm","_lme")
)
head(fit_kimma_all)

# Plot AIC one model vs the other
fit_kimma_all %>% 
  ggplot(aes(x = AIC_lm, y = AIC_lme)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, 
              color="red") +
  coord_fixed()

# Overall, the simple lm is a better fit model. BUT we know the experimental design is paired. You should always consider biology and design first. Then, model fit if those do not clearly define inclusion/exclusion of covariates or random (paired) effects.
# Here, we have a paired design so the paired lme model is preferred.
mean(fit_kimma_all$AIC_lm)
mean(fit_kimma_all$AIC_lme)

# Checkout the end of https://bigslu.github.io/workshops/2024_SEATRAC_series/3_linear_model_rnaseq.html for examples of BIGpicture functions to create more common RNAseq analysis plots.
