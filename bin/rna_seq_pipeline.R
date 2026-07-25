#!/usr/bin/env Rscript

# ==========================================

# RNA-seq DESeq2 Pipeline Script

# ==========================================

suppressPackageStartupMessages({
library(DESeq2)
library(pheatmap)
library(ggplot2)
})

# ==========================================

# 1. Input Arguments

# ==========================================

args <- commandArgs(trailingOnly = TRUE)

count_file   <- args[1]
metadata_file <- args[2]
outdir       <- args[3]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("📂 Count file:", count_file, "\n")
cat("📂 Metadata file:", metadata_file, "\n")

# ==========================================

# 2. Load Count Data

# ==========================================

count <- read.table(count_file,
header = TRUE,
sep = "",
check.names = FALSE)

# Set GeneID as rownames

rownames(count) <- count$GeneID

# Remove GeneID column

count <- count[, -1]

cat("✅ Count matrix loaded:", dim(count), "\n")

# ==========================================

# 3. Load Metadata

# ==========================================

metadata <- read.table(metadata_file,
header = TRUE,
sep = "",
stringsAsFactors = TRUE)

rownames(metadata) <- metadata$sample
metadata <- metadata[, -1]

cat("✅ Metadata loaded:", dim(metadata), "\n")

#---------------------------------------------
cat("\n===== COUNT MATRIX =====\n")
cat("Dimensions:", dim(count), "\n")

cat("\nColumn names first print:\n")
print(colnames(count))

count <- count[, 6:ncol(count)]



# Remove featureCounts suffix
colnames(count) <- sub(
    "\\.fastq\\.Aligned\\.sortedByCoord\\.out\\.bam$",
    "",
    colnames(count)
)
cat("\nColumn names:\n")
print(colnames(count))


cat("\nMetadata:\n")
print(rownames(metadata))

# Check for missing samples
count_samples <- colnames(count)
meta_samples  <- rownames(metadata)

if (!setequal(count_samples, meta_samples)) {

    cat("Samples in count but not metadata:\n")
    print(setdiff(count_samples, meta_samples))

    cat("Samples in metadata but not count:\n")
    print(setdiff(meta_samples, count_samples))

    stop("❌ Sample names do not match")
}

cat("✅ All sample names are present in both files\n")

# Reorder metadata to match count matrix
metadata <- metadata[count_samples, ]

# Final verification
stopifnot(identical(colnames(count), rownames(metadata)))

cat("✅ Metadata reordered to match count matrix\n")

#---------------------------------------------

# ==========================================

# 4. Validation Check

# ==========================================

if (!all(colnames(count) == rownames(metadata))) {
stop("❌ Sample mismatch between count matrix and metadata")
}

cat("✅ Sample names matched\n")

# ==========================================

# 5. Create DESeq2 Dataset

# ==========================================
print("test1")
dds <- DESeqDataSetFromMatrix(
countData = count,
colData = metadata,
design = ~ genotype + zinc
)

# ==========================================

# 6. Filter Low-Count Genes

# ==========================================

dds <- dds[rowSums(counts(dds)) > 10, ]

cat("✅ Genes after filtering:", nrow(dds), "\n")

# ==========================================

# 7. Run DESeq2

# ==========================================

dds <- DESeq(dds)

# ==========================================

# 8. Normalized Counts

# ==========================================

normalized_counts <- counts(dds, normalized = TRUE)
print("test2")
cat("outdir =", outdir, "\n")
print(outdir)
write.csv(normalized_counts,
file = file.path(outdir, "normalized_counts.csv"))

# ==========================================

# 9. Variance Stabilization

# ==========================================

vsd <- varianceStabilizingTransformation(dds)

# ==========================================

# ==========================================
# Create output directories (FIRST)
# ==========================================
dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "qc"), showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)
dir.create(file.path(outdir, "dea"), showWarnings = FALSE)

# ==========================================
# PCA Plot
# ==========================================
pdf(file.path(outdir, "plots/PCA.pdf"))
plotPCA(vsd, intgroup = c("genotype", "zinc"))
dev.off()

# ==========================================
# Sample Distance Heatmap
# ==========================================
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pdf(file.path(outdir, "qc/sample_heatmap.pdf"))
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists)
dev.off()
print("test3")
# ==========================================
# Boxplot
# ==========================================
pdf(file.path(outdir, "qc/boxplot.pdf"))
boxplot(assay(vsd),
        main = "Boxplot of Samples",
        las = 2)
dev.off()
print("test4")
# ==========================================
# DESeq2 Results
# ==========================================
res <- results(dds)

print(resultsNames(dds))
print("test5")
res_genotype <- results(dds, name = "genotype_WT_vs_uzcR")
res_genotype <- res_genotype[order(res_genotype$padj, na.last = TRUE), ]
res_zinc <- results(dds, name = "zinc_Zn_vs_noZn")
res_zinc <- res_zinc[order(res_zinc$padj, na.last = TRUE), ]
print("test6")
dir.create(file.path(outdir, "dea"),
           recursive = TRUE,
           showWarnings = FALSE)

write.csv(
    as.data.frame(res_genotype),
    file = file.path(outdir, "dea", "DE_genotype_results.csv"),
    row.names = TRUE
)
write.csv(
    as.data.frame(res_zinc),
    file = file.path(outdir, "dea", "DE_zinc_results.csv"),
    row.names = TRUE
)
print("test7")
# ==========================================
# Significant Genes
# ==========================================
resSig_genotype <- res_genotype[which(res_genotype$padj < 0.05), ]
resSig_zinc <- res_zinc[which(res_zinc$padj < 0.05), ]

write.csv(as.data.frame(resSig_genotype),
          file = file.path(outdir, "dea/DE_genotype_significant.csv"))
write.csv(as.data.frame(resSig_zinc),
          file = file.path(outdir, "dea/DE_zinc_significant.csv"))
print("test8")
# ==========================================
# Volcano Plot
# ==========================================
pdf(file.path(outdir, "plots/volcano_genotype.pdf"))
plot(res_genotype$log2FoldChange,
     -log10(res_genotype$pvalue),
     pch = 20,
     main = "Volcano Plot",
     xlab = "Log2 Fold Change",
     ylab = "-log10 p-value")
dev.off()
print("test9")
pdf(file.path(outdir, "plots/volcano_zinc.pdf"))
plot(res_zinc$log2FoldChange,
     -log10(res_zinc$pvalue),
     pch = 20,
     main = "Volcano Plot",
     xlab = "Log2 Fold Change",
     ylab = "-log10 p-value")
dev.off()
print("test9")
# ==========================================
# Heatmap of Top Genes
# ==========================================
topGenes_genotype <- head(order(res_genotype$padj), 20)
topGenes_zinc <- head(order(res_zinc$padj), 20)


pdf(file.path(outdir, "plots/top_genes_heatmap_g.pdf"))
pheatmap(assay(vsd)[topGenes_genotype, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE)
dev.off()
pdf(file.path(outdir, "plots/top_genes_heatmap_z.pdf"))
pheatmap(assay(vsd)[topGenes_zinc, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE)
dev.off()
print("test10")
# ==========================================
# DONE
# ==========================================
cat("🎉 RNA-seq pipeline completed successfully!\n")