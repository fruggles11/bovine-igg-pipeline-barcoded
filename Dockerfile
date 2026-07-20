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
    git \
    build-essential \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment with bioinformatics tools
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Activate the micromamba env for all subsequent RUN steps
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install minibar (the calacademy-research ONT/Sanger demultiplexer, not the
# unrelated "minibar" progress-bar package on PyPI)
RUN pip install edlib && \
    git clone https://github.com/calacademy-research/minibar.git /opt/minibar && \
    chmod +x /opt/minibar/minibar.py && \
    ln -s /opt/minibar/minibar.py /usr/local/bin/minibar.py

# Install amplicon_sorter (strip Windows CRLF line endings from the repo's
# source file, which break the #!/usr/bin/env python3 shebang under exec)
RUN pip install amplicon-sorter || \
    (git clone https://github.com/avierstr/amplicon_sorter.git /opt/amplicon_sorter && \
    sed -i 's/\r$//' /opt/amplicon_sorter/amplicon_sorter.py && \
    chmod +x /opt/amplicon_sorter/amplicon_sorter.py && \
    ln -s /opt/amplicon_sorter/amplicon_sorter.py /usr/local/bin/amplicon_sorter.py)

# IgBLAST is installed via environment.yml (bioconda), not a manually
# downloaded binary -- NCBI only ships x64 Linux builds, which have to run
# under Rosetta/QEMU emulation on Apple Silicon and are unreliable there.
# The bioconda build resolves to a native linux-aarch64 binary instead, and
# already bundles internal_data/optional_file under its share directory.
ENV IGDATA=/opt/conda/share/igblast

# Set PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Set working directory
WORKDIR /data

# Default command
CMD ["bash"]
