# ============================================================
# RNA-seq Nextflow Pipeline Docker Image
# Author: Stalin Dany Shaji
# ============================================================

FROM rocker/tidyverse:4.5.1

LABEL maintainer="Stalin Dany Shaji"
LABEL org.opencontainers.image.title="RNA-seq Nextflow Pipeline"
LABEL org.opencontainers.image.description="RNA-seq differential expression analysis using Nextflow, STAR, featureCounts and DESeq2"
LABEL org.opencontainers.image.version="1.0.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:/usr/local/bin:$PATH
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# ============================================================
# System Dependencies
# ============================================================

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    bzip2 \
    ca-certificates \
    build-essential \
    openjdk-17-jdk \
    procps \
    && rm -rf /var/lib/apt/lists/*

# ============================================================
# Install Nextflow
# ============================================================

RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# ============================================================
# Install Miniconda
# ============================================================

RUN wget -qO /tmp/miniconda.sh \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash /tmp/miniconda.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniconda.sh

# ============================================================
# Configure Conda
# ============================================================

RUN conda config --system --add channels conda-forge && \
    conda config --system --add channels bioconda && \
    conda config --system --add channels defaults

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install Python 3.12 and Mamba
RUN conda install -n base -c conda-forge -y \
    python=3.12 \
    mamba && \
    conda clean -afy

# ============================================================
# Install Bioinformatics Tools
# ============================================================

RUN mamba install -y \
    fastqc \
    fastp \
    star \
    subread \
    samtools \
    bedtools \
    gsl \
    pigz \
    parallel && \
    conda clean -afy

# ============================================================
# Install DESeq2
# ============================================================

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install('DESeq2', ask=FALSE, update=FALSE)"

# ============================================================
# Additional R Packages
# ============================================================

RUN R -e "install.packages(c('pheatmap','RColorBrewer'), repos='https://cloud.r-project.org')"

# ============================================================
# Copy Pipeline Scripts
# ============================================================

WORKDIR /workspace

COPY bin/ /workspace/bin/

RUN chmod -R +x /workspace/bin

# ============================================================
# Verify Installation
# ============================================================

RUN java -version

RUN nextflow -version

RUN fastqc --version

RUN fastp --version

RUN STAR --version

RUN featureCounts -v

RUN samtools --version

RUN bedtools --version

RUN Rscript -e "library(DESeq2); sessionInfo()"

# ============================================================
# Default Working Directory
# ============================================================

WORKDIR /workspace

CMD ["/bin/bash"]