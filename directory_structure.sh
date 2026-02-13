#!/bin/bash

# Project root directory
PROJECT="rna-seq-nextflow"

# Create directories
mkdir -p $PROJECT/{configs,modules,bin,docker,data}

# Create main files
touch $PROJECT/main.nf
touch $PROJECT/nextflow.config
touch $PROJECT/README.md

# Create config files
touch $PROJECT/configs/base.config
touch $PROJECT/configs/docker.config
touch $PROJECT/configs/singularity.config
touch $PROJECT/configs/slurm.config

# Create module files
touch $PROJECT/modules/fastqc.nf
touch $PROJECT/modules/trim.nf
touch $PROJECT/modules/align.nf
touch $PROJECT/modules/featurecounts.nf
touch $PROJECT/modules/deseq2.nf
touch $PROJECT/modules/report.nf

# Create bin script
touch $PROJECT/bin/deseq2_analysis.R

# Create docker file
touch $PROJECT/docker/Dockerfile

# Create data file
touch $PROJECT/data/samplesheet.csv

echo "Directory structure for $PROJECT created successfully!"

