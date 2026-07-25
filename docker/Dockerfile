FROM rocker/tidyverse:latest

LABEL maintainer="Stalin Dany Shaji"

# System dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    git \
    openjdk-17-jdk \
    python3 \
    python3-pip \
    build-essential \
    && rm -rf /var/lib/apt/lists/*



# Install Miniconda
ENV CONDA_DIR=/opt/conda

RUN wget -qO /tmp/miniconda.sh \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash /tmp/miniconda.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniconda.sh

ENV PATH=${CONDA_DIR}/bin:$PATH

# Verify conda exists
RUN conda --version

# Configure channels
RUN conda config --system --add channels conda-forge && \
    conda config --system --add channels bioconda && \
    conda config --system --add channels defaults

# Update conda first
#RUN conda update -n base -c defaults conda -y

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install mamba
RUN conda install -n base -c conda-forge mamba -y
RUN mamba install -y \
    -c conda-forge \
    -c bioconda \
    fastqc \
    fastp \
    star \
    subread \
    samtools && \
    conda clean -afy

# Install MultiQC
RUN pip3 install multiqc --break-system-packages

# Install Bioconductor
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c('DESeq2','SummarizedExperiment','BiocParallel'), ask=FALSE, update=FALSE)"

# Install CRAN packages
RUN R -e "install.packages(c('ggplot2','pheatmap','rmarkdown','knitr','RColorBrewer','tidyr','dplyr'), repos='https://cloud.r-project.org')"

WORKDIR /workspace

COPY bin/ bin/

RUN chmod -R +x bin/

# Verify installation
RUN fastqc --version && \
    fastp --version && \
    STAR --version && \
    featureCounts -v && \
    samtools --version && \
    Rscript -e "library(DESeq2)"

CMD ["/bin/bash"]