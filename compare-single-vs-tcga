###########
# Comparing single cell and tcga markers
###########

# Load single cell markers
nc_down <- readRDS("markers/nc_down_markers.rds") 
nc_up <- readRDS("markers/nc_up_markers.rds")

# Load DeSeq tcga markers
tcga_down <- readRDS("markers/tcga_down_markers.rds") 
tcga_up <- readRDS("markers/tcga_up_markers.rds")

# Find common single cell and tcga markers
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

# Keeping only gene names and log2fc values from the data frames
head(nc)
names(nc)[names(nc) == "avg_log2FC"] <- "nc_log2fc"
simple_nc <- nc[c("gene", "nc_log2fc")]
head(simple_nc)

head(tcga)
names(tcga)[names(tcga) == "log2fc"] <- "tcga_log2fc"
simple_tcga <- tcga[c("gene", "tcga_log2fc")]
head(simple_tcga)
print(nrow(simple_tcga))
print(nrow(simple_nc))

# Merge the two data frames based on the gene column
nc_tcga_log2fc <- merge(simple_nc, simple_tcga, by = "gene")
head(nc_tcga_log2fc)
print(nrow(nc_tcga_log2fc))

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

legend_order <- c("Significant in single cell", "Significant in TCGA", "Significant in both", "Not significant")

nc_tcga_log2fc$Label <- nc_tcga_log2fc$gene 
nc_tcga_log2fc$Label[nc_tcga_log2fc$Significance != "Significant in both"] <- NA

ggplot(nc_tcga_log2fc, aes(x = nc_log2fc, y = tcga_log2fc, col=Significance, label=Label)) +
  geom_point(alpha = 0.5) +
  labs(x = "Single cell log2FC", y = "TCGA log2FC", title = "Single cell log2FC vs. TCGA log2FC") +
  scale_color_manual(values=c("#50C878", "blue", "red", "grey"), breaks = legend_order) +
  geom_text_repel() +
  theme_classic()

# Correlation
correlation <- cor(nc_tcga_log2fc$nc_log2fc, nc_tcga_log2fc$tcga_log2fc, use = "complete.obs")
correlation

# making a data frame containing markers significant in both tcga and single cell results
nc_tcga_log2fc$Significance[nc_tcga_log2fc$gene %in% intersect(nc_genes, tcga_genes)] <- "Significant in both"

significant_both <- subset(nc_tcga_log2fc, Significance == "Significant in both")
significant_both

new <- significant_both[, 1:3]
head(new)

# Calculate average log2fc and maximum log2fc so we have two different ways of sorting the data frame
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
