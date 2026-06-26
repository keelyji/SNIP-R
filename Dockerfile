# SNIP-R analysis pipeline
#
# Base image pins R 4.4.1 (the version used to generate the paper figures;
# see notebooks/PartA*.html sessionInfo() block).
#
# Build:    docker build -t snipr:latest .
# Run:      docker run --rm -v "$PWD":/work -w /work snipr:latest \
#               Rscript /opt/snipr/scripts/partC_screen_analysis.R \
#                 --counts example_data/SupplementaryTable4_RawScreenData.xlsx \
#                 --config config/example_partC_config.yaml \
#                 --out-dir output/

FROM rocker/tidyverse:4.4.1

LABEL org.opencontainers.image.title="SNIP-R"
LABEL org.opencontainers.image.description="Reproducible analysis pipeline for SNIP-R CRISPR enhancer screens"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.source="https://github.com/keelyji/SNIP-R"

# ---- system deps (bedtools for Part B FASTA extraction) -------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bedtools \
        python3 \
        python3-pip \
        python3-venv \
        python3-dev \
        && \
    rm -rf /var/lib/apt/lists/*

# ---- R packages ------------------------------------------------------------
# rocker/tidyverse already ships dplyr, tidyr, readr, ggplot2, readxl,
# patchwork, ggrepel. We add the extras used by the CLI scripts and
# bedtoolsr for the FASTA-extraction step in Part B.
RUN R -e "install.packages(c('optparse','yaml','bedtoolsr','ggrepel','patchwork'), repos='https://packagemanager.posit.co/cran/__linux__/jammy/latest')"

# ---- Snakemake (optional orchestration layer) -----------------------------
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --no-cache-dir 'snakemake==7.32.4' 'pulp<2.8' pyyaml
ENV PATH="/opt/venv/bin:${PATH}"

# ---- copy pipeline ---------------------------------------------------------
WORKDIR /opt/snipr
COPY scripts/  /opt/snipr/scripts/
COPY config/   /opt/snipr/config/
COPY Snakefile /opt/snipr/Snakefile

# Convenience wrappers on $PATH so users can call the scripts directly.
RUN ln -s /opt/snipr/scripts/partB_prep_flanks.R     /usr/local/bin/snipr-prep-flanks    && \
    ln -s /opt/snipr/scripts/partB_pair_grnas.R      /usr/local/bin/snipr-pair-grnas     && \
    ln -s /opt/snipr/scripts/partC_screen_analysis.R /usr/local/bin/snipr-screen-analyze && \
    chmod +x /opt/snipr/scripts/*.R

# Default working directory for `docker run -v "$PWD":/work ...`
WORKDIR /work

CMD ["bash"]
