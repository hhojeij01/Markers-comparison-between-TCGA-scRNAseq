# Used the Seurat pipeline (link: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#finding-differentially-expressed-features-cluster-biomarkers)
library(Seurat)

# single cell data comes from BC atlas (link: https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers#study-summary)
# GEO code: GSE176078 
# I logged into the single cell portal (previous BC atlas link) and downloaded data there
# I downloaded matrix, barcodes, and features & uncompressed files via terminal (MacOS) beforehand
# Also downloaded the csv file of the atlas so I could get cell type data

# Load single-cell matrix (saved the matrix as a .rds object beforehand because it took a while to load it otherwise)
matrix <- readRDS("data/matrix.rds")
head(matrix)

# Load in features and barcodes 
barcodes <- read.table("data/barcodes.tsv", header = FALSE)
features <- read.table("data/features.tsv", header = FALSE)
dim(barcodes)
dim(features)

# The matrix we read doesn't have column or row names
# We need to assign barcodes to columns (cells) and features to rows (genes)

# Checking the dimensions of objects
dim(matrix)
dim(features)
dim(barcodes)

# The matrix number of columns match the barcodes
# The matrix number of rows match the features (gene names), BUT we have too many columns in the features frame
# Check those
features

# The features frame has 2 seemingly identical columns, but I'm gonna keep them in because I don't know if they're identical everywhere, since it is a long list
# THe 3rd and 4th column just say "gene" and "expression" for every row, so I'm gonna remove these columns
features <- features[,1:2]

# Now that barcodes (cells) and features (gene names) are ready, we can assign them to the matrix
rownames(matrix) <- features$V1
colnames(matrix) <- barcodes$V1
head(matrix)
dim(matrix)

# Create the seurat object
seurat_object <- CreateSeuratObject(counts = matrix, project = "bc_atlas")
seurat_object

# Checking the name, number, and order of cells 
head(seurat_object)
colnames(seurat_object)
ncells <- Cells(seurat_object)
ncells

# Extracting the major cell types from the study (we need this to be able to analyze cancer and normal cells)
cell_type_annotation <- read.csv("data/Whole_miniatlas_meta.csv")
cell_type_annotation
celltype_major <- cell_type_annotation$celltype_major

# Remove the first row since it doesn't contain any cell information, just tags
celltype_major <- celltype_major[-1]
head(celltype_major)
celltype_major

# Add a column to the metadata in the seurat object called celltype_major
seurat_object <- AddMetaData(seurat_object, metadata = celltype_major, col.name = "celltype_major")
seurat_object@meta.data

# We can check the seurat object's metadata for the first 5 cells
head(seurat_object@meta.data, 5)

# Calculating mitochondrial QC metrics for genes starting with MT
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# We now have 3 quality control metrics stored in metada: nCount_RNA, nFeature_RNA, and percent.mt

# Visualize QC metrics as a violin plot
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1
plot2

# Subset filtering (first to only normal and cancer cells, then get rid of low expression cells)
Idents(seurat_object) <- "celltype_major"
Idents(seurat_object)
seurat_object@meta.data

cells_to_keep <- c("Normal Epithelial", "Cancer Epithelial")
seurat_object <- subset(seurat_object, idents = cells_to_keep)
tail(seurat_object, 200)

seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Check if we only have normal or cancer cells now
head(seurat_object, 200)
tail(seurat_object, 200)

# Unique() gives us the unique "types" we have and colnames() just gives column names 
unique(seurat_object$celltype_major)
colnames(seurat_object@meta.data)
unique(seurat_subset$celltype_major)
head(seurat_subset)
head(seurat_object)

# Normalize 
seurat_object <- NormalizeData(seurat_object)

# Find variable features 
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
head(seurat_object@meta.data, 5)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)
top10

# plot variable features with and without labels
vplot1 <- VariableFeaturePlot(seurat_object)
vplot1
vplot2 <- LabelPoints(plot = vplot1, points = top10, repel = TRUE)
vplot2
vplot1 + vplot2

# Scale the data
# Note: I only scaled the variable genes because using all genes was memory exhaustive
var_genes <- VariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object, features = var_genes)

# PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)

# Determining dimensionality
# Using heatmaps
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
DimPlot(seurat_object, reduction = "pca")
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE)

# Using p values with jackstrawplot
seurat_object <- JackStraw(seurat_object, num.replicate = 100)
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)
JackStrawPlot(seurat_object, dims = 1:15)

# Using elbow plot - based on this, I'll use the first 5 PCs
ElbowPlot(seurat_object)

# Clustering
seurat_object <- FindNeighbors(seurat_object, dims = 1:5)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
head(Idents(seurat_object), 5)

head(cell_type_annotation)

seurat_object$
head(seurat_object)
seurat_object <- RunUMAP(seurat_object, dims = 1:5)
DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident")
DimPlot(seurat_object, reduction = "umap", group.by = "celltype_major")

DimPlot(seurat_object, group.by="seurat_clusters")
seurat_object$seurat_clusters

# Most of previous steps were for clustering & making the umaps

# Save processed seurat_object
saveRDS(seurat_object, file = "data/NC_seurat_object.rds")
test <- readRDS("data/NC_seurat_object.rds")
head(test@meta.data)
head(test)
head(seurat_object)
head(seurat_object@meta.data)
unique(test$celltype_major)

# Finding differentially expressed markers
library(Seurat)
seurat_object <- readRDS("data/NC_seurat_object.rds")
#bc_markers <- FindAllMarkers(seurat_object, logfc.threshold = 0)

# Normal vs. Cancer
markers <- FindMarkers(seurat_object, group.by = "celltype_major", ident.1 = "Normal Epithelial", ident.2 = "Cancer Epithelial", logfc.threshold = 0)
head(markers)

# CAFs vs. Cancer
# markers <- FindMarkers(seurat_object, group.by = "celltype_major", ident.1 = "CAFs", ident.2 = "Cancer Epithelial", logfc.threshold = 0)
# head(markers)

# Normal and CAFs vs. Cancer
# markers <- FindMarkers(seurat_object, group.by = "celltype_major", ident.1 = c("Normal Epithelial", "CAFs"), ident.2 = "Cancer Epithelial", logfc.threshold = 0)
# head(markers)

# Save markers as rds object
saveRDS(markers, file = "data/NC_markers.rds")

# Plotting the log2fc against pval 
plot(markers$avg_log2FC, -log10(markers$p_val))
plot(de$avg_log2FC, -log10(de$p_val))
head(de)

# Making a dataframe from the markers
de <- readRDS(file = "data/NC_markers.rds")

# Remove markers with p-vals of 0
de <- de[de$p_val != 0, ]

# Create a new column called diffexpressed and assign "NO" to all markers (i.e. not differentially expressed, which we'll change in a sec)
de$diffexpressed <- "NO"
head(de)

# Plot to see the range
plot(de$avg_log2FC, -log10(de$p_val))
range(de$p_val)

# if log2Foldchange > 2.5 and pvalue < 0.05, set as "UP" for upregulated
de$diffexpressed[de$avg_log2FC > 0.5 & de$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -2.5 and pvalue < 0.05, set as "DOWN" for downregulated
de$diffexpressed[de$avg_log2FC < -0.5 & de$p_val_adj < 0.05] <- "DOWN"

# Add column for labels (gonna use this for ggplot labels, this column is just gonna be the gene names but only present for diffexpressed)
de$delabel <- NA

# Add column for genes (genes are the rownames but it was a pain filtering using rownames)
de$gene <- rownames(de)
head(de)

# Only add labels (gene names) for the differentially expressed markers, in other words if diffexpressed = 0 there won't be a label in ggplot
de$delabel[de$diffexpressed != "NO"] <- de$gene[de$diffexpressed != "NO"]
head(de)
tail(de)

# plot with colours & labels using ggplot2
library(ggplot2)
library(ggrepel)
ggplot(de, aes(x = avg_log2FC, y = -log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 1) +
  theme_bw() +
  labs(x = "log2FC", y = "-log10(p_val_adj)", title = "Volcano Plot: P-value vs. log2FC of markers") +
  geom_vline(xintercept=c(-0.5, 0.5), col="red", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "longdash") +
  scale_x_continuous(breaks = seq(-5, 5, 1)) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() +
  ylim(0,400)

# Save the data frame for plotting
saveRDS(de, file = "data/NC_complete_markers.rds")
test <- readRDS("data/NC_complete_markers.rds")
head(test)
plot(test$avg_log2FC, -log10(test$p_val))
plot(test$avg_log2FC, -log10(test$p_val_adj))
head(test)

# Dataset stats
unique(seurat_object$celltype_major)
normal_counts <- table(seurat_object$celltype_major)["Normal Epithelial"]
cancer_counts <- table(seurat_object$celltype_major)["Cancer Epithelial"]
cafs <- table(seurat_object$celltype_major)["CAFs"]

# Print the number of normal cells
print(normal_counts)
print(cancer_counts)
print(cafs)

nrow(seurat_object)
tail(seurat_object)

normal_counts <- table(sub$celltype_major)["Normal Epithelial"]
cancer_counts <- table(sub$celltype_major)["Cancer Epithelial"]
unique(sub$celltype_major)

# Check the original upregulated and downregulated counts
markers <- readRDS("data/NC_complete_markers.rds")
head(markers)
dim(markers)
plot(markers$avg_log2FC, -log10(markers$p_val_adj))

down_count <- table(markers$diffexpressed)["DOWN"]
print(down_count)

up_count <- table(markers$diffexpressed)["UP"]
print(up_count)

# Extract downregulated list of markers
nc_down <- markers[markers$diffexpressed == "DOWN", ]
nc_down_markers <- rownames(nc_down)

# Extract upregulated list of markers
nc_up <- markers[markers$diffexpressed == "UP", ]
nc_up_markers <- rownames(nc_up)

# Save them
saveRDS(nc_down_markers, file="markers/nc_down_markers.rds")
saveRDS(nc_up_markers, file="markers/nc_up_markers.rds")

# Make sure it's all good
down <- readRDS("markers/nc_down_markers.rds") 
up <- readRDS("markers/nc_up_markers.rds") 
down
up
