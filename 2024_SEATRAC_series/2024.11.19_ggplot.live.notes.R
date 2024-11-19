#Setup
## Rproject
## Rscript
## data (dat_voom.RData)

#### Libraries ####
library(tidyverse)
library(limma)

#### Data ####
load("data/dat_voom.RData")

#### ggplot layers ####
# data
# aesthetics (aes)
# geometries (geom)
# statistics (stat)
# guides (annotation)
# facets

#### simplest plot ####
# data + geom + aes
ggplot(data = dat$targets) +
  aes(y = total_seq) +
  geom_boxplot()

# Add x variable
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot()

#### Multiple geoms ####
# Add individual dots
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot() +
  geom_point()

# jitter points to avoid overlap
## use a seed for reproducibility
set.seed(42)
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2)

# order of layers matters
## dots under boxplot
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition)+
  geom_jitter(width=0.2) +
  geom_boxplot(outlier.shape = NA)

## dots on top of boxplot
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.2)

# Add horizontal cutoff line
## Again think about order and where you want the line
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition)+
  geom_jitter(width=0.2) +
  geom_hline(yintercept = 1E7)+
  geom_boxplot(outlier.shape = NA)

#### Setting aesthetics ####
# aes priority
# space 
# color
# shape/size
# other things

# Add color
ggplot(data = dat$targets) +
  aes(y = total_seq, x = condition, color=condition)+
  geom_jitter(width=0.2) +
  geom_hline(yintercept = 1E7)+
  geom_boxplot(outlier.shape = NA)

#### PCA ####
#Calculate PCA
dat.pca <- prcomp(t(dat$E), scale. = TRUE,
                  center = TRUE)
dat.pca$x[1:3,1:3]

#Format PCA data
dat.pca.meta <- as.data.frame(dat.pca$x) %>% 
  rownames_to_column("libID") %>% 
  inner_join(dat$targets, by="libID")

# plot pca
## Use color
ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, color=condition)+
  geom_point(size=5)

## Use shape
ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, shape=condition)+
  geom_point(size=5)

## Use color and shape
ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, color=condition,
      shape=condition)+
  geom_point(size=5)

#### Changing the look (theme) ####
# get rid of grey (preset theme)
ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, color=condition,
      shape=condition)+
  geom_point(size=5) +
  theme_classic()

ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, color=condition,
      shape=condition)+
  geom_point(size=5) +
  theme_bw()

# get rid of grey (one item in theme)
ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, color=condition,
      shape=condition)+
  geom_point(size=5) +
  theme(panel.background = element_rect(color="black", fill = "white"))

#### Change labels ####
# Get PCA percent explained
summary(dat.pca)$importance

ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, color=condition,
      shape=condition)+
  geom_point(size=5) +
  theme_bw() +
  labs(x = "PC1 (24%)", y = "PC2 (11%)")
# ylab() xlab() also exist

#### Coordinate systems ####
# force x and y ratio
ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, color=condition,
      shape=condition)+
  geom_point(size=5) +
  theme_bw() +
  labs(x = "PC1 (24%)", y = "PC2 (11%)",
       color = "Mtb infection",
       shape = "Mtb infection") +
  coord_fixed(ratio = 1)

#### Change legend labels ####
# by changing the data
dat.pca.meta %>% 
  mutate(cond2 = case_when(
    condition=="Media"~"Control",
    TRUE~condition)) %>% 
  ggplot() +
  aes(x=PC1, y=PC2, color=cond2,
      shape=cond2)+
  geom_point(size=5) +
  theme_bw() +
  labs(x = "PC1 (24%)", y = "PC2 (11%)",
       color = "Mtb infection",
       shape = "Mtb infection") +
  coord_fixed(ratio = 1)

# by changing the scale
ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, color=condition)+
  geom_point(size=5) +
  theme_bw() +
  labs(x = "PC1 (24%)", y = "PC2 (11%)",
       color = "Mtb infection") +
  coord_fixed(ratio = 1) +
  scale_color_discrete(labels =c("Control","Mtb"))

#### move the legend #####
ggplot(dat.pca.meta) +
  aes(x=PC1, y=PC2, color=condition)+
  geom_point(size=5) +
  theme_bw() +
  labs(x = "PC1 (24%)", y = "PC2 (11%)",
       color = "Mtb infection") +
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom",
        legend.direction = "vertical")

#### volcano plot ####
#Load the data
load("data/kimma_result.RData")
fit_kimma$lme %>% View

ggplot(fit_kimma$lme) +
  aes(x=estimate, y=-log10(FDR)) +
  geom_point()

#### facet ####
# Multiple plots of subsets of data
ggplot(fit_kimma$lme) +
  aes(x=estimate, y=-log10(FDR)) +
  geom_point() +
  facet_wrap(~variable, scale="free")

#### Another aes example ####
#add and change color
plot1 <- fit_kimma$lme %>% 
  filter(variable == "condition") %>% 
  #create color group for significant up and down
  mutate(col.group = case_when(
    FDR < 0.01 & estimate > 0 ~"up",
    FDR < 0.01 & estimate < 0 ~ "down",
    TRUE~"NS"
  )) %>% 
  mutate(col.group = factor(col.group, 
                            levels=c("up","down","NS"))) %>% 
  ggplot() +
  aes(x=estimate, y=-log10(FDR),
      color=col.group) +
  geom_point() +
  scale_color_manual(values = c("down"="blue",
                                "NS"="grey",
                                "up"="red")) +
  geom_hline(yintercept = -log10(0.01),
             lty="dashed") +
  lims(x=c(-9,9)) +
  theme_bw(base_size = 20) +
  labs(color="Fold change") 
plot1

#### Change axis ticks and lines ####
plot1 +
  scale_x_continuous(breaks = c(seq(-10,10, 1)), limits = c(-9,9))

#### Label points ####
library(ggrepel)

# Too many!
plot1 +
  geom_text_repel(aes(label=gene))


plot2 <- fit_kimma$lme %>% 
  filter(variable == "condition") %>% 
  mutate(col.group = case_when(
    FDR < 0.01 & estimate > 0 ~"up",
    FDR < 0.01 & estimate < 0 ~ "down",
    TRUE~"NS"
  )) %>% 
  mutate(col.group = factor(col.group, 
                            levels=c("up","down","NS"))) %>% 
  mutate(lab = case_when(estimate>8 ~ gene,
                         TRUE~ NA)) %>% 
  ggplot() +
  aes(x=estimate, y=-log10(FDR),
      color=col.group) +
  geom_point() +
  scale_color_manual(
    values = c("down"="blue",
               "NS"="grey",
               "up"="red")) +
  geom_hline(yintercept = -log10(0.01),
             lty="dashed") +
  lims(x=c(-9,9)) +
  theme_bw(base_size = 20) +
  labs(color="Fold change") +
  geom_text_repel(aes(label=lab), 
                  max.overlaps = Inf,
                  show.legend = FALSE)

#### Save a plot ####
ggsave(plot2, file="volcano.png", width=8, height=5)

#### Resources for choosing colors ####
# https://davidmathlogic.com/colorblind
# https://colorbrewer2.org/
# viridis
# https://github.com/thomasp85/scico
# https://coolors.co/

library(viridis)
plot2 +
  aes(color=estimate) +
  scale_color_viridis(option = "magma")

#### My cheatsheet #####
#Not just ggplot. https://kdillmcfarland.github.io/r/R-tips/

#### add greek symbol pdf ####
#on Mac
#install cairo (I was wrong. Not Xcode, it's homebrew)
# https://formulae.brew.sh/formula/cairo

# Use that package to save
ggsave(plot2, file="volcano.pdf", width=8, height=5, device = cairo_pdf)

# Or use expression instead of symbol code
labs(y = expression("X ("*Delta*"Ct)"))
