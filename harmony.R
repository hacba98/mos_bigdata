# This method works very similar to normal intergration pipeline
# Only need to modified data object with provided Harmony function

# Read and combine 2 Mouse retina raw dataset in 'mouse' folder
sub1.data <- readRDS('mouse/sub1.rds')
sub1 <- as.Seurat(sub1.data, counts = "counts")
rm(sub1.data)

sub2.data <- readRDS('mouse/sub2.rds')
sub2 <- as.Seurat(sub2.data, counts = "counts")
rm(sub2.data)

gc()

raw.combined <- merge(sub1, y=sub2, add.cell.ids=c("Sub1","Sub2"), project="Mouse Retina")

# the Harmony algorithm iteratively corrects PCA embeddings.
raw.combined <- FindVariableFeatures(raw.combined)
raw.combined <- Seurat::ScaleData(raw.combined, features = GetVariableFeatures(raw.combined, arg), verbose = FALSE)
raw.combined <- RunPCA(raw.combined, arg, logger, seed)
# pass the Seurat object and specify which variable(s) to integrate out
raw.combined <- harmony::RunHarmony(raw.combined, "orig.ident", plot_convergence = FALSE, verbose = FALSE, epsilon.harmony = -Inf, max.iter.harmony = 30)
DimPlot(raw.combined, reduction = "umap", label = TRUE, pt.size = .1)