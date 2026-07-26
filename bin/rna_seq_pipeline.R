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
design <- as.formula(args[3])
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
if(!"sample" %in% colnames(metadata)){
    stop("Metadata should contain a column named 'sample'.")
}
rownames(metadata) <- metadata$sample
metadata <- metadata[, colnames(metadata) != "sample"]

cat("✅ Metadata loaded:", dim(metadata), "\n")

cat("Metadata columns:\n")
print(colnames(metadata))

cat("Design formula:\n")
print(design)

# Extract variable names from the design formula
design_vars <- all.vars(design)

cat("Variables in design:\n")
print(design_vars)

# Find missing variables
missing_vars <- setdiff(design_vars, colnames(metadata))

if (length(missing_vars) > 0) {
    stop(
        paste0(
            "\n❌ Invalid design formula.\n",
            "Missing metadata column(s): ",
            paste(missing_vars, collapse = ", "),
            "\n\nAvailable metadata columns:\n",
            paste(colnames(metadata), collapse = ", ")
        )
    )
}

cat("✅ Design formula validated successfully.\n")

#---------------------------------------------
cat("\n===== COUNT MATRIX =====\n")
cat("Dimensions:", dim(count), "\n")

# Set gene IDs as row names
rownames(count) <- count$Geneid

cat("\nColumn names first print:\n")
print(colnames(count))

sample_cols <- !colnames(count) %in% c(
    "Geneid","GeneID",
    "Chr","Start","End","Strand","Length"
)

count <- count[, sample_cols]

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
design = design
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

coef_names <- resultsNames(dds)

print(coef_names)

dir.create(
    file.path(outdir, "dea"),
    recursive = TRUE,
    showWarnings = FALSE
)

for (coef in coef_names) {

    # Skip intercept
    if (coef == "Intercept")
        next

    res <- results(dds, name = coef)

    res <- res[order(res$padj, na.last = TRUE), ]

    write.csv(
        as.data.frame(res),
        file = file.path(outdir, "dea", paste0(coef, ".csv")),
        row.names = TRUE
    )

    # ------------------------------
    # Significant genes
    # ------------------------------
    resSig <- res[which(res$padj < 0.05), ]

    write.csv(
        as.data.frame(resSig),
        file = file.path(outdir,
                         "dea",
                         paste0(coef, "_significant.csv")),
        row.names = TRUE
    )

    # ------------------------------
    # Volcano Plot
    # ------------------------------
    pdf(file.path(outdir,
                  "plots",
                  paste0("volcano_", coef, ".pdf")))

    plot(
        res$log2FoldChange,
        -log10(res$pvalue),
        pch = 20,
        main = coef,
        xlab = "Log2 Fold Change",
        ylab = "-log10(p-value)"
    )

    dev.off()

    # ------------------------------
    # Top 20 genes
    # ------------------------------
    topGenes <- head(order(res$padj), 20)

    pdf(file.path(outdir,
                  "plots",
                  paste0("heatmap_", coef, ".pdf")))

    pheatmap(
        assay(vsd)[topGenes, ],
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        main = coef
    )

    dev.off()
}

cat("🎉 RNA-seq pipeline completed successfully!\n")


