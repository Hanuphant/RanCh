# Use an older miniconda base image to better support legacy packages
FROM continuumio/miniconda3:4.12.0

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Install system dependencies needed for R and bioinformatics tools
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libgsl-dev \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy conda environment file
COPY RanCh/betaenvironment.yml .

# Create conda environment with special handling for older packages
# Note: Using mamba for faster dependency resolution with older packages
RUN conda install mamba -n base -c conda-forge && \
    mamba env create -f betaenvironment.yml && \
    conda clean -afy

# Activate environment and install additional R packages if needed
RUN echo "source activate beta" > ~/.bashrc
ENV PATH /opt/conda/envs/beta/bin:$PATH
SHELL ["conda", "run", "-n", "beta", "/bin/bash", "-c"]

# Install additional R packages that might be missing
RUN conda run -n beta R -e "if (!require('shiny', quietly = TRUE)) install.packages('shiny', repos='https://cran.rstudio.com/')"
RUN conda run -n beta R -e "if (!require('shinydashboard', quietly = TRUE)) install.packages('shinydashboard', repos='https://cran.rstudio.com/')"
RUN conda run -n beta R -e "if (!require('DT', quietly = TRUE)) install.packages('DT', repos='https://cran.rstudio.com/')"
RUN conda run -n beta R -e "if (!require('plotly', quietly = TRUE)) install.packages('plotly', repos='https://cran.rstudio.com/')"

# Install BiocManager and DESeq2 if not already available
RUN conda run -n beta R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cran.rstudio.com/')"
RUN conda run -n beta R -e "if (!require('DESeq2', quietly = TRUE)) BiocManager::install('DESeq2')"

# Create necessary directories
RUN mkdir -p /app/www /app/Analysis/BETA_Demo

# Set permissions for BETA tool (from cistrome_beta package)
RUN chmod +x /opt/conda/envs/beta/bin/BETA 2>/dev/null || echo "BETA binary not found or already executable"

# Expose port for Shiny app
EXPOSE 3838

# Create startup script for interactive use
RUN echo '#!/bin/bash\nsource activate beta\nexec "$@"' > /app/entrypoint.sh && \
    chmod +x /app/entrypoint.sh

# Set the entrypoint to activate conda environment
ENTRYPOINT ["/app/entrypoint.sh"]

# Default command opens an interactive bash shell with conda environment activated
CMD ["/bin/bash"]

