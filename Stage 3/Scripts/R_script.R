#setwd
setwd('/Users/HP/Desktop/RScripts')

#Load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)


# Load data
counts <- read.delim("counts2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
meta <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = TRUE)
meta$condition <- factor(meta$Phase, levels = c("Chronic Infection", "Acute Infection"))

# 2. Prepare metadata 
# Make sure the metadata has a column named 'sample'
# and set it as rownames (so DESeq2 can match samples)
rownames(meta) <- meta$Sample

# 3. Keep only count columns 
# Drop annotation columns like Chr, Start, End, Strand, Length
counts_mat <- counts[, !(colnames(counts) %in% c("Chr", "Start", "End", "Strand", "Length"))]

# 4. Add Geneid as rownames 
rownames(counts_mat) <- counts_mat$Geneid
counts_mat <- counts_mat[, -1]  # remove Geneid column

# 5. Match order of samples
counts_mat <- counts_mat[, rownames(meta)]

# 6. Sanity check 
cat("\n✅ Samples in counts:", length(colnames(counts_mat)))
cat("\n✅ Samples in metadata:", length(rownames(meta)), "\n")

if (!all(colnames(counts_mat) == rownames(meta))) {
  stop("❌ Sample names in counts and metadata do not match. Please check!")
}

# 7. Preview 
cat("\nPreview (samples only):\n")
print(head(counts_mat[, 1:5]))

cat("\nPreview (with Geneid as rownames):\n")
print(head(counts_mat))


# 8. Create DESeq2 dataset 
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = meta,
                              design = ~ Phase) 

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

#preview
dds
dds$Sample
dds$Phase

#Perform the differential expression analysis
dds <- DESeq(dds)

# 4. Results
res <- results(dds, contrast = c("Phase", "Acute Infection", "Chronic Infection"))
res <- res[order(res$padj), ]


#look at your result
head(res)
summary(res)
# Save full results table
write.csv(as.data.frame(res), file = "DESeq2_sa.csv", row.names = TRUE)


# rlog for visualization
rld <- rlog(dds, blind = FALSE)

# 5. PCA
png("PCA_sa.png", width = 6, height = 5, units = "in", res = 300)
plotPCA(rld, intgroup = "Phase") + ggtitle("PCA: Acute vs Chronic")
dev.off()

# 6. Sample distance heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- rownames(meta)
png("SampleDistanceHeatmap_sa.png", width = 6, height = 6, units = "in", res = 300)
pheatmap(sampleDistMatrix, annotation_col = meta)
dev.off()

# 7. Volcano plot (ggplot2)
volcano_data <- as.data.frame(res)
volcano_data$gene <- rownames(volcano_data)

volcano_data$threshold <- "Not Sig"
volcano_data$threshold[volcano_data$padj < 0.05 & volcano_data$log2FoldChange > 1] <- "Up"
volcano_data$threshold[volcano_data$padj < 0.05 & volcano_data$log2FoldChange < -1] <- "Down"

png("Volcano_sa.png", width = 6, height = 5, units = "in", res = 300)
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.6, size = 1.5, na.rm = TRUE) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  theme_minimal() +
  xlab("log2 Fold Change (Acute vs Chronic)") +
  ylab("-log10 Adjusted p-value") +
  ggtitle("Volcano Plot: Acute vs Chronic PJI")
dev.off()

# 8. Heatmap of top 30 DEGs
topGenes <- head(order(res$padj), 30)
png("Heatmap_sa_top30.png", width = 7, height = 8, units = "in", res = 300)
pheatmap(assay(rld)[topGenes, , drop = FALSE],
         annotation_col = meta,
         show_rownames = TRUE, cluster_cols = TRUE)
dev.off()


sig <- subset(res, !is.na(padj) & padj < 0.05)
up <- subset(sig, log2FoldChange > 1)
down <- subset(sig, log2FoldChange < -1)
write.csv(up, file = "upregulated.csv")
write.csv(down, file ="downregulated.csv")

# 9. Export normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts), "NormalizedCounts_sa.csv", row.names = TRUE)
