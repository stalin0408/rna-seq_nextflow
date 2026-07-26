# RNA-Seq Differential Expression Pipeline

A reproducible RNA-Seq analysis pipeline built with **Nextflow DSL2**. The pipeline performs quality control, read trimming, alignment, gene quantification, and differential expression analysis using DESeq2. It is designed to be portable through Docker and configurable for different experimental designs.

---

## Workflow

FASTQ
→ FastQC
→ FastP
→ FastQC
→ STAR Alignment
→ FeatureCounts
→ DESeq2
→ QC Plots & Differential Expression Results

---

## Prerequisites

- Nextflow (v25 or later)
- Docker
- Reference genome
- Gene annotation (GTF)

---

## Reference Genome

Download the following files from **Ensembl**.

### Genome FASTA

```
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

### Gene Annotation

```
Homo_sapiens.GRCh38.115.gtf
```

Place them inside:

```
genome/
├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
└── Homo_sapiens.GRCh38.115.gtf
```

---

## Input Data

Place all FASTQ files inside the `data/` directory.

Example dataset:

```
SRR4280548.fastq.gz
SRR4280549.fastq.gz
SRR4280550.fastq.gz
SRR4280551.fastq.gz
SRR4280552.fastq.gz
SRR4280553.fastq.gz
SRR4280554.fastq.gz
SRR4280555.fastq.gz
SRR4280556.fastq.gz
```

The included `metadata.tsv` corresponds to this example dataset.

---

## Metadata Format

The metadata file must contain a sample identifier and one or more experimental variables.

Example:

| sample | genotype | zinc |
|---------|----------|------|
| SRR4280548 | WT | Zn |
| SRR4280549 | WT | noZn |
| SRR4280550 | uzcR | Zn |
| ... | ... | ... |

The **sample** column must match the FASTQ filenames (without the `.fastq.gz` extension).

---

## Running the Pipeline

Basic execution:

```bash
nextflow run main.nf \
    -resume \
    -with-docker
```

Specify any valid DESeq2 design formula:

```bash
nextflow run main.nf \
    -resume \
    -with-docker \
    --design "~ genotype + zinc"
```

Another example:

```bash
nextflow run main.nf \
    --design "~ condition + batch"
```

The variables used in `--design` **must exist as column names in `metadata.tsv`**.

---

## Output

```
results/
├── fastqc/
├── fastp/
├── star/
├── featurecounts/
└── deseq2_results/
    ├── normalized_counts.csv
    ├── dea/
    ├── plots/
    └── qc/
```

The pipeline automatically generates:

- FastQC reports
- FastP reports
- STAR alignment outputs
- Gene count matrix
- Normalized counts
- Differential expression tables
- Significant gene lists
- PCA plot
- Volcano plots
- Sample distance heatmap
- Top differentially expressed gene heatmaps

---

## Customisable Parameters

| Parameter | Description |
|-----------|-------------|
| `--design` | DESeq2 design formula |
| `--input_dir` | Input FASTQ directory |
| `--meta` | Metadata file |
| `--gtf` | Annotation file |
| `--genome` | STAR genome index |
| `--min_quality` | FASTP quality threshold |
| `--min_length` | Minimum read length |

---

## Docker

Build the Docker image:

```bash
docker build -t rnaseq-pipeline:1.0 -f docker/Dockerfile .
```

Run with Docker:

```bash
nextflow run main.nf \
    -profile docker \
    --design "~ genotype + zinc"
```

---

## Citation

If you use this pipeline, please cite the original software packages:

- Nextflow
- FastQC
- FastP
- STAR
- FeatureCounts (Subread)
- DESeq2

---