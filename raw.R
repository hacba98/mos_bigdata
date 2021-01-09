# Require library Seurat.
# Require memory > 32 Gb RAM
# Can be installed by command: 'install.packages("Seurat")''
# Or follow Seurat document
library('Seurat')

# Read and combine 2 Mouse retina raw dataset in 'mouse' folder
sub1.data <- readRDS('mouse/sub1.rds')
sub1 <- as.Seurat(sub1.data, counts = "counts")
rm(sub1.data)

sub2.data <- readRDS('mouse/sub2.rds')
sub2 <- as.Seurat(sub2.data, counts = "counts")
rm(sub2.data)

gc()

raw.combined <- merge(sub1, y=sub2, add.cell.ids=c("Sub1","Sub2"), project="Mouse Retina")

# Data preprocessing
raw.combined <- NormalizeData(raw.combined, normalization.method = "LogNormalize", scale.factor = 10000)
raw.combined <- FindVariableFeatures(raw.combined)

# Run dimensional reduction - PCA
raw.combined <- RunPCA(raw.combined, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

# Graph-based clustering
raw.combined <- FindNeighbors(raw.combined, reduction = "pca", dims = 1:75, nn.eps = 0.5)
raw.combined <- FindClusters(raw.combined, resolution = 3, n.start = 10)

# Visualize with tsne & UMAP
raw.combined <- RunTSNE(raw.combined, dims = 1:75, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)
raw.combined <- RunUMAP(raw.combined, dims = 1:75, min.dist = 0.75)

library(ggplot2)
p1 <- DimPlot(raw.combined, reduction = "tsne", pt.size = 0.1) + ggtitle(label = "FIt-SNE")
p2 <- DimPlot(raw.combined, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP")
p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)
(p1 + p2) & NoLegend()

browser()