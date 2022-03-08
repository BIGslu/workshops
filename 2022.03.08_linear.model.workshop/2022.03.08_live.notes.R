#### Packages ####
# Install packages
install.packages(c("tidyverse", "lme4", "car", "BiocManager","devtools"))
BiocManager::install("limma")
devtools::install_github("BIGslu/kimma")

# Load packages
library(tidyverse)
library(broom)
library(lme4)
library(car)
library(limma)
library(kimma)

#### Data ####
# Load and format data into 1 data frame with 1 gene's expression
dat <- as.data.frame(example.voom$E) %>% 
  rownames_to_column("geneName") %>% 
  pivot_longer(-geneName, names_to = "libID") %>% 
  inner_join(example.voom$targets, by = "libID") %>% 
  inner_join(example.voom$genes, by = "geneName") %>% 
  filter(hgnc_symbol == "ZNF439")

dat

#### Syntax linear model in R ####
# y = x in "normal" writing is y ~ x in R

# Formula are a specific class of object in R
class(y ~ x)
class(y = x) #error

#### T-test ####
?t.test #help in RStudio

# Fit model and estimate P-values
TT <- t.test(formula = value ~ virus, data = dat)

# Extract t-test result in nice data frame
TT.df <- broom::tidy(TT)
TT.df$p.value

#### ANOVA ####
# Fit model and estimate P-values
AV <- aov(formula = value ~ virus, data = dat)
# Extract results
AV.df <- broom::tidy(AV)

# P-values will be slightly different due to P-value estimation algorithm
TT.df$p.value
AV.df$p.value[1]

# Fit another model which would not have worked with a t-test
t.test(formula = value ~ virus + asthma, data = dat) #error
tidy(aov(formula = value ~ virus + asthma, data = dat))

#### Linear model ####
# Fit model and estimate P-values
LM <- lm(formula = value ~ virus + asthma, data = dat)
# Extract results
tidy(LM)

# Remove intercept in a model
# Be careful. Just because you have only categorical variables like here,
# doesn't mean you should remove the intercept
LM <- lm(formula = value ~ 0 + virus + asthma, data = dat)
tidy(LM)

#### More on reference levels ####
# R assumes you reference level is the first alphabetically
sort(unique(dat$asthma))

# unless you code a factor with specific levels
unique(dat$virus)

#### Assumptions ####
# Balance? Yes
dat %>% count(virus)
dat %>% count(asthma)

# Variance equal? Yes
dat %>%  group_by(virus) %>%  summarise(variance = var(value))

# Independent? No
dat %>% count(donorID)

#### Linear mixed effects model ####
# fixed and random effects
# lme4 package
# Other popular option nlme

# Fit model
LME <- lmer(formula = value ~ virus + asthma + (1|donorID), data = dat)

# Estimate P-value and extract results
tidy(LME) #error
tidy(car::Anova(LME))

#### Scaling up ####
# Packages limma, dream,kimma
# Also DESeq2 but not shown here as the data formatting is very different
# and it does not support random effects

## Limma
# Make model matrix. This is limma's version of a formula
mm <- model.matrix(~ virus + asthma, data = example.voom$targets)
# Fit linear regression. Similar to lm( )
fit <- lmFit(object = example.voom$E, design = mm)
# Estimate p-values. Similar to tidy( ) or Anova( )
efit <- eBayes(fit)
# Get results in a data frame
fdr <- extract_lmFit(design = mm, fit = efit)

# 3 rows for 3 variables (intercept, virus, asthma) for each of 1000 genes
nrow(fdr)

##dream R package
# See full workshop notes

##kimma
fit3 <- kmFit(dat = example.voom, patientID = "donorID",
              model = "~ virus + asthma + (1|donorID)", 
              run.lme = TRUE, run.lm = TRUE)

fdr3 <- fit3$lme

#### Fin ####