# Require library Seurat.
# Require memory > 32 Gb RAM
# Can be installed by command: 'install.packages("Seurat")'
# Or follow Seurat document

# Read and combine 2 Mouse retina raw dataset in 'mouse' folder
sub1.data <- readRDS('mouse/sub1.rds')
sub1 <- as.Seurat(sub1.data, counts = "counts")
rm(sub1.data)

sub2.data <- readRDS('mouse/sub2.rds')
sub2 <- as.Seurat(sub2.data, counts = "counts")
rm(sub2.data)

gc()

# Prior to finding anchors, we perform standard preprocessing (log-normalization),
# and identify variable features individually for each.

pancreas.list <- c(sub1, sub2)
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst",
        nfeatures = 2000, verbose = FALSE)
}

reference.list <- pancreas.list[c("Sub1", "Sub2")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

# visualization with new integrated matrix
# run PCA and visualize result with UMAP
library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE,
    repel = TRUE) + NoLegend()
p1 + p2
