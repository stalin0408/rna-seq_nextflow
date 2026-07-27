# RNA-Seq Differential Expression Pipeline

A reproducible **RNA-Seq differential expression analysis pipeline** built using **Nextflow DSL2**. The pipeline performs quality assessment, read trimming, genome alignment, gene quantification, differential expression analysis, and consolidated quality reporting through **MultiQC**.

The workflow is containerized using Docker to ensure portability and reproducibility across different computing environments.

---

## Workflow

```
FASTQ
   │
   ▼
FastQC
   │
   ▼
FastP
   │
   ▼
STAR Alignment
   │
   ▼
FeatureCounts
   │
   ▼
DESeq2
   │
   ▼
MultiQC
```

---

## Features

- Nextflow DSL2 modular workflow
- Docker support for reproducible execution
- FastQC quality assessment
- FastP read trimming and filtering
- STAR genome alignment
- Gene quantification using FeatureCounts
- Differential expression analysis using DESeq2
- Configurable experimental design formula
- Consolidated quality reports using MultiQC
- Publication-ready visualizations

---

## Prerequisites

- Nextflow (v25 or later)
- Docker

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

The provided `metadata.tsv` corresponds to this example dataset.

---

## Metadata Format

The metadata file must contain a **sample** column and one or more experimental variables.

Example:

| sample | genotype | zinc |
|---------|----------|------|
| SRR4280548 | WT | Zn |
| SRR4280549 | WT | noZn |
| SRR4280550 | uzcR | Zn |
| ... | ... | ... |

The **sample** column must match the FASTQ filenames **without the `.fastq.gz` extension**.

---

## Running the Pipeline

### Basic Execution

```bash
nextflow run main.nf \
    -resume \
    -with-docker
```

### Specify a DESeq2 Design Formula

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

The variables used in the design formula **must exist as column names in `metadata.tsv`**.

---

## Output Structure

```
results/
├── fastqc/
├── fastp/
├── star/
├── featurecounts/
├── deseq2_results/
│   ├── normalized_counts.csv
│   ├── dea/
│   ├── plots/
│   └── qc/
└── multiqc/
    └── multiqc_report.html
```

---

## Generated Reports

### Quality Control

- FastQC reports
- FastP trimming reports
- STAR alignment statistics
- FeatureCounts assignment summary
- MultiQC consolidated report

### Differential Expression

- Raw gene count matrix
- Normalized count matrix
- Differential expression results
- Significant gene lists
- PCA plot
- Volcano plots
- Sample distance heatmap
- Top differentially expressed gene heatmaps

---

## Configurable Parameters

| Parameter | Description |
|-----------|-------------|
| `--design` | DESeq2 design formula |
| `--input_dir` | Input FASTQ directory |
| `--meta` | Metadata file |
| `--gtf` | Gene annotation file |
| `--genome` | STAR genome index |
| `--min_quality` | FASTP quality threshold |
| `--min_length` | Minimum read length |

---

## Docker

Build the Docker image:

```bash
docker build -t rnaseq-pipeline:1.0 -f docker/Dockerfile .
```

Run the pipeline:

```bash
nextflow run main.nf \
    -profile docker \
    --design "~ genotype + zinc"
```

---

## Software

| Tool | Version |
|------|---------|
| Nextflow | DSL2 |
| FastQC | Latest supported |
| FastP | Latest supported |
| STAR | Latest supported |
| FeatureCounts (Subread) | Latest supported |
| DESeq2 | Latest supported |
| MultiQC | v1.35 |

---

## Citation

If you use this pipeline in your research, please cite the original publications of:

- Nextflow
- FastQC
- FastP
- STAR
- FeatureCounts (Subread)
- DESeq2
- MultiQC

---

## Author

**Stalin Dany Shaji**

Bioinformatics Engineer