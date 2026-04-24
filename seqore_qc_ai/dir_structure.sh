#!/bin/bash

PROJECT="qc_filter_pipeline"

echo "Creating project structure..."

mkdir -p $PROJECT/{bin,modules/qc,modules/filtering,modules/post_qc,workflows,data,results/qc,results/filtered,results/post_qc,logs}

cd $PROJECT || exit

# Core files
touch main.nf nextflow.config run_pipeline.sh

# Bin scripts
touch bin/detect_read_type.py
touch bin/parse_fastqc.py
touch bin/parse_nanoplot.py

# QC modules
touch modules/qc/fastqc.nf
touch modules/qc/nanoplot.nf

# Filtering modules
touch modules/filtering/fastp.nf
touch modules/filtering/filtlong.nf

# Post-QC modules
touch modules/post_qc/fastqc_post.nf
touch modules/post_qc/nanoplot_post.nf

# Workflows
touch workflows/long_read.nf
touch workflows/short_read.nf

echo "Structure created successfully!"