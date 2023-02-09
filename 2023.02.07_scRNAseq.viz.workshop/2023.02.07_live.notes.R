#Packages
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(circlize)

#Load data
load("pbmc_clean.RData")

#Violin plot
##Percent mitochondrial seqs

pbmc[["percent.mt"]]
class(pbmc[["percent.mt"]])

##Seurat version
VlnPlot(pbmc, features = "percent.mt", 
        group.by = "orig.ident")

##ggplot version
ggplot(data = pbmc[["percent.mt"]]) +
  aes(x = 1, y = percent.mt) +
  geom_point()

ggplot(data = pbmc[["percent.mt"]]) +
  aes(x = 1, y = percent.mt) +
  geom_jitter(height = 0, width = 0.2) +
  geom_violin()

ggplot(data = pbmc[["percent.mt"]]) +
  aes(x = 1, y = percent.mt) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.4)

##final ggplot
set.seed(765)
p1 <- ggplot(data = pbmc[["percent.mt"]]) +
  aes(x = 1, y = percent.mt) +
  geom_violin(fill = "salmon") + #Add color to violin
  geom_jitter(height = 0, width = 0.4, size = 0.5) + #Make points smaller
  #Add nice labels
  labs(x = "Identity pbmc3k", 
       y = "",
       title = "percent.mt") +
  #Format background and axes
  theme_classic() +
  #Remove meaningless x-axis values
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

## Save
ggsave(filename = "mito.pct.violin.png", 
       plot = p1, width = 2, height = 3)


#Ordination
##UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

### Seurat version
DimPlot(pbmc, reduction = "umap", 
        group.by = "orig.ident")

DimPlot(pbmc, reduction = "umap", 
        group.by = "ident")

### ggplot version
umap_xy <- as.data.frame(pbmc[["umap"]]@cell.embeddings)
head(umap_xy)

cell_clusters <- as.data.frame(pbmc$seurat_clusters)
head(cell_clusters)

umap_xy_clusters <- umap_xy %>% 
  rownames_to_column("cell") %>% 
  full_join(
    cell_clusters %>% 
      rownames_to_column("cell")
  ) %>% 
  rename(ident = `pbmc$seurat_clusters`)

head(umap_xy_clusters)

### plot!
ggplot(data = umap_xy_clusters) +
  aes(x = UMAP_1, y = UMAP_2) +
  geom_point(size = 0.5) +
  aes(color = ident) +
  theme_classic()

## Gene expression UMAP
#### Seurat version
FeaturePlot(pbmc, reduction = "umap",
            features = c("LYZ","MS4A1"))

### ggplot version
#Define gene of interest
genes_of_interest <- c("LYZ", "MS4A1")

#Extract normalized gene expression data
pbmc_counts <- as.data.frame(GetAssayData(
  object = pbmc, slot = "data"))
pbmc_counts[1:3,1:3]

#Filter for gene of interest
umap_xy_expression <- pbmc_counts %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% genes_of_interest) %>% 
  #long format
  pivot_longer(-gene, names_to = "cell") %>% 
  #combine with UMAP data
  full_join(umap_xy_clusters)
head(umap_xy_expression)

ggplot(data = umap_xy_expression) +
  aes(x = UMAP_1, y = UMAP_2) +
  geom_point(size = 0.5) +
  aes(color = value) + #Color by expression
  scale_color_gradient(low = "lightgrey", 
                       high = "blue") + #Change colors
  facet_wrap(~ gene, scales = "free") + #Split by gene
  labs(color = "") +
  theme_classic()

#Heatmap
head(pbmc_markers)
top3 <- pbmc_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_log2FC)
head(top3)

## Seurat version
DoHeatmap(pbmc, features = top3$gene)

DoHeatmap(pbmc, features = top3$gene,
          disp.max = 10)

## ggplot version
pbmc_counts_log <- as.data.frame(GetAssayData(
  pbmc, slot = "scale.data"))

#Order cells by clusters. Remember that we originally 
# got these data from as.data.frame(pbmc$seurat_clusters)
cell_order <- cell_clusters %>% 
  arrange(`pbmc$seurat_clusters`)

pbmc_counts_log_top3 <- pbmc_counts_log %>% 
  #Filter top genes
  rownames_to_column("gene") %>% 
  filter(gene %in% top3$gene) %>% 
  #Remove outliers
  mutate_if(is.numeric, ~ifelse(. > 2.5, 2.5, .)) %>% 
  #Reorder cells
  select(gene, all_of(rownames(cell_order))) %>% 
  #Reorder genes
  mutate(gene = factor(gene, levels = top3$gene)) %>% 
  arrange(gene) %>% 
  #Convert to matrix
  column_to_rownames("gene") %>% 
  as.matrix()

pbmc_counts_log_top3[1:2,1:2]

Heatmap(pbmc_counts_log_top3,
        use_raster = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE)

Heatmap(pbmc_counts_log_top3,
        use_raster = FALSE,
        cluster_rows = TRUE,
        cluster_columns = TRUE)

