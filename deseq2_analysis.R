library(DESeq2)

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]
meta_file <- args[2]



counts <- read.delim(counts_file, comment.char="#")
rownames(counts) <- counts$Geneid
counts <- counts[, 7:ncol(counts)]
colnames(counts) = gsub(".fastq.Aligned.sortedByCoord.out.bam", "", colnames(counts))


# Read metadata
coldata <- read.delim(meta_file, row.names=1)

# Ensure sample order matches
counts <- counts[, rownames(coldata)]

# Create DESeq object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ genotype + zinc + genotype:zinc
)

dds <- DESeq(dds)

res <- results(dds)

write.csv(as.data.frame(res), "deseq2_results.csv")

