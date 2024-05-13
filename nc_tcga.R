library(DESeq2)

# Load data
data1 <- readRDS("/Users/hadihojeij/Documents/bc_atlas/DeSeqData/TCGA-BRCA-eset.rds")
data1

# Checking what our data contains
rowData(data1)$gene_id
rowData(data1)$gene_name
colnames(rowData(data1))

names(assays(data1))
colData(data1)

# Check samples
table(colData(data1)$sample_type)

# Keep only tumour and normal samples
data <- data1[, colData(data1)$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
unique(data@colData$sample_type)

# Remove non protein coding genes (rows)
data <- data[rowData(data)$gene_type == "protein_coding", ]
unique(rowData(data)$gene_type)

# Remove gene duplicates
data <- data[!duplicated(rowData(data)$gene_name), ]

# Replace ensembl ID rownames as gene names
rownames(data) <- rowData(data)$gene_name
rowData(data)

# Create counts matrix from unstranded
studyExpr <- assays(data)[["unstranded"]]
head(studyExpr[,1:3])

# Remove genes with NA values or 0 variance (using 1 in function indicates rows/genes)
rvar <- apply(studyExpr, 1, var)
studyExpr <- studyExpr[complete.cases(rvar), ]
rvar <- apply(studyExpr, 1, var)
studyExpr <- studyExpr[rvar != 0, ]

# Remove patients with NA values or 0 variance (using 2 in function indicates columns/patients)
cvar <- apply(studyExpr, 2, var)
studyExpr <- studyExpr[, complete.cases(cvar)]
cvar <- apply(studyExpr, 2, var)
studyExpr <- studyExpr[, cvar != 0]

# Create clinical metadata for colData
clinical <- colData(data)
colnames(clinical)
unique(clinical$sample_type)

# Create Deseq data set from matrix and collect results
dds <- DESeqDataSetFromMatrix(countData = studyExpr, colData = clinical, design = ~sample_type)
dds <- estimateSizeFactors(dds)
#studyExpr <- counts(dds, normalized = TRUE)
test <- DESeq(dds)
results <- results(test)
head(results)
plot(results$log2FoldChange, -log10(results$padj))

# Create data frame
resultsDF <- data.frame(
  gene = rownames(results),
  log2fc = results$log2FoldChange,
  pvalue = results$pvalue
)
head(resultsDF)

# Plot
resultsDF$diffexpressed <- "NO"
resultsDF

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
resultsDF$diffexpressed[resultsDF$log2fc > 2.5 & resultsDF$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
resultsDF$diffexpressed[resultsDF$log2fc < -2.5 & resultsDF$pvalue < 0.05] <- "DOWN"
resultsDF$delabel <- NA
resultsDF$delabel[resultsDF$diffexpressed != "NO"] <- resultsDF$gene[resultsDF$diffexpressed != "NO"]

library(ggplot2)
library(ggrepel)
ggplot(resultsDF, aes(x = log2fc, y = -log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 1) +
  theme_bw() +
  labs(x = "log2FC", y = "-log10(p_val)", title = "Volcano Plot: P-value vs. log2FC of all markers") +
  geom_vline(xintercept=c(-2.5, 2.5), col="red", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "longdash") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() #+
#xlim(-10,10)

# Save data frame
head(resultsDF)
saveRDS(resultsDF, file = "DeSeqData/tcga_results.rds")
head(resultsDF)
test <- readRDS("DeSeqData/tcga_results.rds")
head(test)
plot(test$log2fc, -log10(test$pvalue))


# Check the original upregulated and downregulated counts
markers <- readRDS("DeSeqData/tcga_results.rds")
head(markers)
plot(markers$log2fc, -log10(markers$pvalue))

down_count <- table(markers$diffexpressed)["DOWN"]
print(down_count)

up_count <- table(markers$diffexpressed)["UP"]
print(up_count)

# Extract downregulated
tcga_down <- markers[markers$diffexpressed == "DOWN", ]
tcga_down_markers <- tcga_down$gene
tcga_down_markers

# Extract upregulated
tcga_up <- markers[markers$diffexpressed == "UP", ]
tcga_up_markers <- tcga_up$gene
tcga_up_markers

# Save them
saveRDS(tcga_down_markers, file="markers/tcga_down_markers.rds")
saveRDS(tcga_up_markers, file="markers/tcga_up_markers.rds")

# Make sure it's all good
down <- readRDS("markers/tcga_down_markers.rds") 
up <- readRDS("markers/tcga_up_markers.rds")
down
up

###########
# Compare
###########

# Load single cell markers
nc_down <- readRDS("markers/nc_down_markers.rds") 
nc_up <- readRDS("markers/nc_up_markers.rds")

# Load DeSeq tcga markers
tcga_down <- readRDS("markers/tcga_down_markers.rds") 
tcga_up <- readRDS("markers/tcga_up_markers.rds")

# Single cell and tcga markers
common_down_nc_tcga <- intersect(nc_down, tcga_down)
common_down_nc_tcga

common_up_nc_tcga <- intersect(nc_up, tcga_up)
common_up_nc_tcga

####################################################
# Log2fc scatterplots
#####################################################

# Load in differentially expressed markers data frames
nc <- readRDS("/Users/hadihojeij/Documents/bc_atlas/data/NC_complete_markers.rds")
tcga <- readRDS("/Users/hadihojeij/Documents/bc_atlas/DeSeqData/tcga_results.rds")

head(nc)
names(nc)[names(nc) == "avg_log2FC"] <- "nc_log2fc"
simple_nc <- nc[c("gene", "nc_log2fc")]
head(simple_nc)

head(tcga)
names(tcga)[names(tcga) == "log2fc"] <- "tcga_log2fc"
simple_tcga <- tcga[c("gene", "tcga_log2fc")]
head(simple_tcga)

# Merge the two data frames based on the gene column
nc_tcga_log2fc <- merge(simple_nc, simple_tcga, by = "gene")
head(nc_tcga_log2fc)

# Lists of significant markers and combine down and up
nc_down <- readRDS("markers/nc_down_markers.rds") 
nc_up <- readRDS("markers/nc_up_markers.rds")
nc_genes <- c(nc_down, nc_up)
nc_genes

tcga_down <- readRDS("markers/tcga_down_markers.rds") 
tcga_up <- readRDS("markers/tcga_up_markers.rds") 
tcga_genes <- c(tcga_up, tcga_down)
tcga_genes

# Create the third column with initial value as "not significant"
nc_tcga_log2fc$Significance <- "Not significant"
head(nc_tcga_log2fc)

# Label genes based on significance
nc_tcga_log2fc$Significance[nc_tcga_log2fc$gene %in% nc_genes] <- "Significant in single cell"
nc_tcga_log2fc$Significance[nc_tcga_log2fc$gene %in% tcga_genes] <- "Significant in TCGA"
nc_tcga_log2fc$Significance[nc_tcga_log2fc$gene %in% intersect(nc_genes, tcga_genes)] <- "Significant in both"
head(nc_tcga_log2fc)
nc_tcga_log2fc
unique(nc_tcga_log2fc$Significance)

# Plot
library(ggplot2)
library(ggrepel)
ggplot(nc_tcga_log2fc, aes(x = nc_log2fc, y = tcga_log2fc, col=Significance, label=Label)) +
  geom_point(alpha = 0.5) +
  labs(x = "Single cell log2FC", y = "TCGA log2FC", title = "Single cell log2FC vs. TCGA log2FC") +
  scale_color_manual(values=c("#50C878", "blue", "red", "grey"), breaks = legend_order) +
  geom_text_repel() +
  theme_classic()

legend_order <- c("Significant in single cell", "Significant in TCGA", "Significant in both", "Not significant")

nc_tcga_log2fc$Label <- nc_tcga_log2fc$gene 
nc_tcga_log2fc$Label[nc_tcga_log2fc$Significance != "Significant in both"] <- NA

# Correlation
correlation <- cor(nc_tcga_log2fc$nc_log2fc, nc_tcga_log2fc$tcga_log2fc, use = "complete.obs")
correlation

#
nc_tcga_log2fc$Significance[nc_tcga_log2fc$gene %in% intersect(nc_genes, tcga_genes)] <- "Significant in both"

significant_both <- subset(nc_tcga_log2fc, Significance == "Significant in both")
significant_both

new <- significant_both[, 1:3]
head(new)

# Calculate average log2fc and maximum log2fc
new$nc_log2fc
new$tcga_log2fc
new$avg_log2FC <- rowMeans(new[, c("nc_log2fc", "tcga_log2fc")])

new$max_log2FC <- apply(new[, c("nc_log2fc", "tcga_log2fc")], 1, max)
head(new)

# Sort by avg log2fc
sorted_markers <- new[order(new$avg_log2FC, decreasing = TRUE), ]
head(sorted_markers)

# Sort by max log2fc
sorted_markers <- new[order(new$max_log2FC, decreasing = TRUE), ]
sorted_markers

