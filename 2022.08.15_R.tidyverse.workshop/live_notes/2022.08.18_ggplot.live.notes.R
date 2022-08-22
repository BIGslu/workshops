## Some code has been condensed to reduce redundancy. Suggest running iteratively
## to follow online video i.e. run line 1, then 1+2, then 1+2+3, etc within a ggplot 
## connected by +

# Load packages
library(tidyverse)

# Load data
load("data/dat_voom.RData")

# Build a minimal ggplot

ggplot(data = dat$targets) + # Make a ggplot object (grey box) and link the data you want to use
  aes(y = total_seq) +       # Map aesthetics (aes) to tell what part of the data to plot
  geom_boxplot()             # Give a geom for how the data should be represented

## Split the data along the x axis with an addtl aesthetic
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot()

## aesthetics can be nested within other functions like ggplot() or the geom
ggplot(data = dat$targets,
       aes(y = total_seq, x = condition)) +
  geom_boxplot()

ggplot(data = dat$targets) +
  geom_boxplot(aes(y = total_seq, x = condition))

## Change variable ordering with a factor
### Note how we pipe %>% the modified data into the ggplot and do not have to use data =  inside ggplot()
dat$targets %>% 
  mutate(condition_fct = factor(condition,
                                levels = c("Mtb","Media"))) %>% 
  ggplot() +
  aes(y = total_seq, x = condition_fct) +
  geom_boxplot()
  
## Add dots for each data point on top of boxplot
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_point() + 
  geom_boxplot()

## layer order matters. Here the boxplot covers up some data points
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot() + 
  geom_point()

## jitter so you can see all the data points
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot() + 
  geom_jitter(height = 0, width = 0.2)

## Oh no! There are duplicate data points for the boxplot outlier samples and the jitter points
## Remove boxplot outliers
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.2)

## force jitter location with a seed
set.seed(345786)

ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.2)

set.seed(345786)

ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.2)

# RNAseq plots
## PCA
## Calculate PCA
dat.pca <- prcomp(dat$E)
## Want values per library not per gene so transpose expression data with t()
dat.pca <- prcomp(t(dat$E), scale.=TRUE, center=TRUE)

## PC values are in x
class(dat.pca$x)
dat.pca$x

## Extract PC values and merge with library metadata
dat.pca.meta <- as.data.frame(dat.pca$x) %>% 
  rownames_to_column("libID") %>% 
  inner_join(dat$targets)

## plot PCA
ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2) +
  geom_point()

## color
ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point()
## Again you can put aes nested instead
ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2) +
  geom_point(aes(color = condition))
  
## shape
ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, shape = condition) +
  geom_point()

## labels
### Get % variation explained by each PC
summary(dat.pca)

plot2 <- ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection")

### You can also individually call each label as an argument
ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  xlab("PC1 (24.4%)") +
  ylab("PC2 (11.3%)") +
  ggtitle("PCA of gene expression")

## coordinate system to force spacing on x and y to be the same
## Doesn't change this plot much but can is 1 PC has a different range than the other
ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed()

## themes
ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed() +
  theme_classic()

ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed() +
  theme_bw()

ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed() +
  theme_minimal()

ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed() +
  theme_dark()
##ggthemes package for more themes!!

## add axis at 0,0 
## This breaks some of the tenants of ggplot grammar so it is difficult to do
## I would do something like this, leaving the x and y labels on the edges
ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed() +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## But you can really force it by drawing the new axes with ticks and labels
## Note recommended because the data (points) are mixed in with axis labels (numbers) and it can be confusing

## Create data frame with tick line info
axis_data <- data.frame(
  x = c(-100, -50, 0, 50, 100,
        0, 0, 0, 0, 0),
  xend = c(-100, -50, 0, 50, 100,
           1, 1, 1, 1, 1),
  y = c(0, 0, 0, 0, 0,
        -100, -50, 0, 50, 100),
  yend = c(1, 1, 1, 1, 1,
           -100, -50, 0, 50, 100)
)

ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed() +
  ## All of this just for the axes change
  theme_void() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(angle = 90)) +
  geom_segment(data = axis_data, ## Note how we call a new data object for just this layer
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "black") +
  ## x-axis
  geom_text(data = filter(axis_data, x != 0),
            aes(x = xend, y = yend+2, label = x),
            color = "black") +
  ## y-axis
  geom_text(data = filter(axis_data, y != 0),
            aes(x = xend+2, y = yend, label = y),
            color = "black")

## Paired points by patient

ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  geom_line(color = "black") +
  aes(group = ptID) + ## Note adding group so ggplot knows which points to connect
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed() +
  theme_classic()

## Paired by condition ellipses
ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed() +
  stat_ellipse(level = 0.9) # confidence interval ellipses

## Volcano plot
## Linear model results
load("data/kimma_result.RData")

## Linear model was
## ~ condition + sex + (1|ptID)
head(fit$lme)

## Basic volcano plot
fit$lme %>% 
  filter(variable == "condition") %>% 
  ggplot() +
  aes(x = estimate, y = -log10(FDR)) +
  geom_point()

## Color
fit$lme %>% 
  filter(variable == "condition") %>% 
  ggplot() +
  aes(x = estimate, y = -log10(FDR),
      color = -log10(FDR)) +
  geom_point()

## Better color by significance cutoff
fit$lme %>% 
  filter(variable == "condition") %>% 
  mutate(significance = case_when(
    FDR < 0.01 ~ "FDR < 0.01",
    TRUE ~ "NS")) %>% 
  ggplot() +
  aes(x = estimate, y = -log10(FDR),
      color = significance) +
  geom_point() +
  geom_hline(yintercept = -log10(0.01))

## Addtl significance group for up vs down regulated
fitC <- fit$lme %>% 
  filter(variable == "condition") %>% 
  mutate(significance = case_when(
    FDR < 0.01 & estimate < 0 ~ "down",
    FDR < 0.01 & estimate > 0 ~ "up",
    TRUE ~ "NS"))

ggplot(fitC) +
  aes(x = estimate, y = -log10(FDR),
      color = significance) +
  geom_point() +
  geom_hline(yintercept = -log10(0.01))

## Change colors
ggplot(fitC) +
  aes(x = estimate, y = -log10(FDR),
      color = significance) +
  geom_point() +
  geom_hline(yintercept = -log10(0.01)) +
  scale_color_manual(values = c("blue","grey","red"))

## Scale axis
ggplot(fitC) +
  aes(x = estimate, y = -log10(FDR),
      color = significance) +
  geom_point() +
  geom_hline(yintercept = -log10(0.01)) +
  scale_color_manual(values = c("blue","grey","red")) +
  lims(x = c(-9,9))

## label points (genes)
library(ggrepel)

## Get top 2 most signif up and down genes
top_deg <- fit$lme %>% 
  filter(variable == "condition") %>% 
  mutate(direction = case_when(
    estimate < 0 ~ "down",
    estimate > 0 ~ "up")) %>% 
  group_by(direction) %>% 
  slice_min(FDR, n = 2) %>% 
  pull(gene) ## Gets just the 1 column as a vector
 
top_deg

## Or list specific genes
## Commented out to not run at this time
# top_deg <- c("IFNG","CIITA")

## To be clear, we go back to the original fit data instead of fitC
fit$lme %>% 
  filter(variable == "condition") %>% 
  mutate(significance = case_when(
    FDR < 0.01 & estimate < 0 ~ "down",
    FDR < 0.01 & estimate > 0 ~ "up",
    TRUE ~ "NS")) %>% 
  mutate(label = case_when(gene %in% top_deg ~ gene)) %>% 

ggplot() +
  aes(x = estimate, y = -log10(FDR),
      color = significance) +
  geom_point() +
  geom_hline(yintercept = -log10(0.01)) +
  scale_color_manual(values = c("blue","grey","red")) +
  lims(x = c(-9,9)) +
  geom_text_repel() +
  aes(label = label)
## Note the warning about missing data. Here was want this because we do not want labels for all of the genes

## Beautify
fit$lme %>% 
  filter(variable == "condition") %>% 
  mutate(significance = case_when(
    FDR < 0.01 & estimate < 0 ~ "down",
    FDR < 0.01 & estimate > 0 ~ "up",
    TRUE ~ "NS")) %>% 
  mutate(label = case_when(gene %in% top_deg ~ gene))  %>% 
  
  ggplot() +
  aes(x = estimate, y = -log10(FDR)) +
  geom_point(aes(color = significance)) + # move color to within this geom so it only applies to points, not the later text
  geom_hline(yintercept = -log10(0.01)) +
  scale_color_manual(values = c("blue","grey","red")) +
  lims(x = c(-9,9)) +
  geom_text_repel(show.legend = FALSE, # Remove lettering legend since it's redundant with points
                  fontface = "italic") + # put labels in italics
  aes(label = label) +
  theme_bw()

## facets 
plot1 <- fit$lme %>% 
  filter(variable %in% c("condition", "sex")) %>% # now keep two variables
  mutate(significance = case_when(
    FDR < 0.01 & estimate < 0 ~ "down",
    FDR < 0.01 & estimate > 0 ~ "up",
    TRUE ~ "NS")) %>% 
  mutate(label = case_when(gene %in% top_deg ~ gene))  %>% 
  
  ggplot() +
  aes(x = estimate, y = -log10(FDR)) +
  geom_point(aes(color = significance)) +
  geom_hline(yintercept = -log10(0.01)) +
  scale_color_manual(values = c("blue","grey","red")) +
  geom_text_repel(show.legend = FALSE) +
  aes(label=label) +
  facet_wrap(~ variable, scales = "free") #facet
plot1

## save plots
ggsave(filename = "volcano.plot.png", plot=plot1,
       width=6, height=4)

## Multiple panels not in facet
install.packages("patchwork")
library(patchwork)

### Save two plots
plot2 <- ggplot(dat.pca.meta) +
  aes(x = PC1, y = PC2, color = condition) +
  geom_point() +
  labs(x = "PC1 (24.4%)", y = "PC2 (11.3%)",
       color = "Infection") +
  coord_fixed() +
  theme_classic()
  
plot1 + plot2 + plot_layout(widths=c(2,1))

# Reuse your own theme
my_theme <- theme(...)
my_plot + my_theme
