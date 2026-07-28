#!/usr/bin/env Rscript

# ============================================================
# RNA-seq DESeq2 Pipeline Script
# ============================================================

suppressPackageStartupMessages({

    library(DESeq2)
    library(pheatmap)
    library(ggplot2)

})

# ============================================================
# 1. Read Input Arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

count_file    <- args[1]
metadata_file <- args[2]
design        <- as.formula(args[3])
outdir        <- args[4]

# Create output directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("📂 Count file:", count_file, "\n")
cat("📂 Metadata file:", metadata_file, "\n")

# ============================================================
# 2. Load Count Matrix
# ============================================================

count <- read.table(
    count_file,
    header = TRUE,
    sep = "",
    check.names = FALSE
)

# Set GeneID as row names
rownames(count) <- count$GeneID

# Remove GeneID column
count <- count[, -1]

cat("✅ Count matrix loaded:", dim(count), "\n")

# ============================================================
# 3. Load Sample Metadata
# ============================================================

metadata <- read.table(
    metadata_file,
    header = TRUE,
    sep = "",
    stringsAsFactors = TRUE
)

# Ensure metadata contains a sample column
if (!"sample" %in% colnames(metadata)) {

    stop("Metadata should contain a column named 'sample'.")

}

# Set sample names as row names
rownames(metadata) <- metadata$sample

# Remove sample column
metadata <- metadata[, colnames(metadata) != "sample"]

cat("✅ Metadata loaded:", dim(metadata), "\n")

cat("Metadata columns:\n")

cat("Design formula:\n")

# ============================================================
# 4. Validate Design Formula
# ============================================================

# Extract variable names from the design formula
design_vars <- all.vars(design)

cat("Variables in design:\n")

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

# ============================================================
# 5. Prepare Count Matrix
# ============================================================

cat("\n===== COUNT MATRIX =====\n")
cat("Dimensions:", dim(count), "\n")

# Set gene IDs as row names
rownames(count) <- count$Geneid

cat("\nColumn names (before processing):\n")

# Select only sample columns
sample_cols <- !colnames(count) %in% c(
    "Geneid",
    "GeneID",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length"
)

count <- count[, sample_cols]

# Remove featureCounts suffix from sample names
colnames(count) <- sub(
    "\\.fastq\\.Aligned\\.sortedByCoord\\.out\\.bam$",
    "",
    colnames(count)
)

cat("\nColumn names (after processing):\n")

cat("\nMetadata sample names:\n")

# ============================================================
# 6. Validate Sample Names
# ============================================================

count_samples <- colnames(count)
meta_samples  <- rownames(metadata)

if (!setequal(count_samples, meta_samples)) {

    cat("Samples in count but not metadata:\n")

    cat("Samples in metadata but not count:\n")

    stop("❌ Sample names do not match")

}

cat("✅ All sample names are present in both files\n")

# Reorder metadata to match count matrix
metadata <- metadata[count_samples, ]

# Final verification
stopifnot(identical(colnames(count), rownames(metadata)))

cat("✅ Metadata reordered to match count matrix\n")

# ============================================================
# 7. Final Validation
# ============================================================

if (!all(colnames(count) == rownames(metadata))) {

    stop("❌ Sample mismatch between count matrix and metadata")

}

cat("✅ Sample names matched\n")

# ============================================================
# 8. Create DESeq2 Dataset
# ============================================================

dds <- DESeqDataSetFromMatrix(
    countData = count,
    colData   = metadata,
    design    = design
)

# ============================================================
# 9. Filter Low-Count Genes
# ============================================================

dds <- dds[rowSums(counts(dds)) > 10, ]

cat("✅ Genes after filtering:", nrow(dds), "\n")

# ============================================================
# 10. Run DESeq2 Analysis
# ============================================================

dds <- DESeq(dds)

# ============================================================
# 11. Export Normalized Counts
# ============================================================

normalized_counts <- counts(dds, normalized = TRUE)

cat("outdir =", outdir, "\n")

write.csv(
    normalized_counts,
    file = file.path(outdir, "normalized_counts.csv")
)

# ============================================================
# 12. Variance Stabilizing Transformation (VST)
# ============================================================

vsd <- varianceStabilizingTransformation(dds)

# ============================================================
# 13. Create Output Directories
# ============================================================

dir.create(
    outdir,
    showWarnings = FALSE
)

dir.create(
    file.path(outdir, "qc"),
    showWarnings = FALSE
)

dir.create(
    file.path(outdir, "plots"),
    showWarnings = FALSE
)

dir.create(
    file.path(outdir, "dea"),
    showWarnings = FALSE
)

# ============================================================
# 14. Principal Component Analysis (PCA)
# ============================================================

pdf(file.path(outdir, "plots/PCA.pdf"))

plotPCA(
    vsd,
    intgroup = c("genotype", "zinc")
)

dev.off()

# ============================================================
# 15. Sample Distance Heatmap
# ============================================================

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)

pdf(file.path(outdir, "qc/sample_heatmap.pdf"))

pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists
)

dev.off()

# ============================================================
# 16. Sample Distribution Boxplot
# ============================================================

pdf(file.path(outdir, "qc/boxplot.pdf"))

boxplot(
    assay(vsd),
    main = "Boxplot of Samples",
    las  = 2
)

dev.off()
# ============================================================
# 17. Differential Expression Analysis
# ============================================================

coef_names <- resultsNames(dds)

dir.create(
    file.path(outdir, "dea"),
    recursive = TRUE,
    showWarnings = FALSE
)

for (coef in coef_names) {

    # Skip the intercept coefficient
    if (coef == "Intercept")
        next

    # Extract differential expression results
    res <- results(
        dds,
        name = coef
    )

    # Sort by adjusted p-value
    res <- res[order(res$padj, na.last = TRUE), ]

    # ========================================================
    # Export complete results
    # ========================================================

    write.csv(
        as.data.frame(res),
        file = file.path(
            outdir,
            "dea",
            paste0(coef, ".csv")
        ),
        row.names = TRUE
    )

    # ========================================================
    # Export significantly differentially expressed genes
    # ========================================================

    resSig <- res[which(res$padj < 0.05), ]

    write.csv(
        as.data.frame(resSig),
        file = file.path(
            outdir,
            "dea",
            paste0(coef, "_significant.csv")
        ),
        row.names = TRUE
    )

    # ========================================================
    # Volcano Plot
    # ========================================================

    pdf(
        file.path(
            outdir,
            "plots",
            paste0("volcano_", coef, ".pdf")
        )
    )

    plot(
        res$log2FoldChange,
        -log10(res$pvalue),
        pch  = 20,
        main = coef,
        xlab = "Log2 Fold Change",
        ylab = "-log10(p-value)"
    )

    dev.off()

    # ========================================================
    # Heatmap of Top 20 Differentially Expressed Genes
    # ========================================================

    topGenes <- head(order(res$padj), 20)

    pdf(
        file.path(
            outdir,
            "plots",
            paste0("heatmap_", coef, ".pdf")
        )
    )

    pheatmap(
        assay(vsd)[topGenes, ],
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        main = coef
    )

    dev.off()

}

# ============================================================
# 18. Pipeline Completed Successfully
# ============================================================

cat("🎉 RNA-seq pipeline completed successfully!\n")