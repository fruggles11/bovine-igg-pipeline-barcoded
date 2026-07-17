# Dockerfile for Bovine IgG Repertoire Analysis Pipeline
FROM mambaorg/micromamba:1.5.1

LABEL maintainer="fruggles11"
LABEL description="Bovine IgG repertoire analysis pipeline for Oxford Nanopore data"

# Set up micromamba environment
USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment with bioinformatics tools
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Activate the micromamba env for all subsequent RUN steps
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install minibar
RUN pip install minibar

# Install amplicon_sorter
RUN pip install amplicon-sorter || \
    (git clone https://github.com/avierstr/amplicon_sorter.git /opt/amplicon_sorter && \
    chmod +x /opt/amplicon_sorter/amplicon_sorter.py && \
    ln -s /opt/amplicon_sorter/amplicon_sorter.py /usr/local/bin/amplicon_sorter.py)

# Set up IgBLAST
RUN mkdir -p /opt/igblast && \
    cd /opt/igblast && \
    wget -q https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-linux.tar.gz && \
    tar -xzf ncbi-igblast-1.22.0-x64-linux.tar.gz && \
    rm ncbi-igblast-1.22.0-x64-linux.tar.gz && \
    ln -s /opt/igblast/ncbi-igblast-1.22.0/bin/* /usr/local/bin/

# Copy IgBLAST internal data
RUN cp -r /opt/igblast/ncbi-igblast-1.22.0/internal_data /usr/local/share/igblast/ || true
RUN cp -r /opt/igblast/ncbi-igblast-1.22.0/optional_file /usr/local/share/igblast/ || true

ENV IGDATA=/usr/local/share/igblast

# Set PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Set working directory
WORKDIR /data

# Default command
CMD ["bash"]
