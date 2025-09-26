# 1. Start from a trusted, official R base image
FROM rocker/r-ver:4.3.2

# 2. Install ALL required system-level dependencies.
#    - libglpk-dev: The GNU Linear Programming Kit, required by the 'igraph' R package.
#    - libfftw3-dev: For Fast Fourier Transform, required by Rnmr1D dependencies.
#    - git, libxml2-dev, etc.: Other dependencies we've already identified.
RUN apt-get update -qq && apt-get install -y \
    libarchive-dev  \
    libglpk-dev \
    libfftw3-dev \
    git \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. Install all required R packages from CRAN and Bioconductor in one go.
#    This is more efficient for Docker builds.
RUN R -e '                                                               \
    install.packages(c(                                                  \
        "archive",                                                       \
        "BiocManager", "devtools", "plumber", "drc", "jsonlite",         \
        "ggplot2", "base64enc"                                           \
    ));                                                                      \
    BiocManager::install(c(                                              \
        "impute", "MassSpecWavelet", "pcaMethods"                        \
    ), update=FALSE);                                                    \
'

# 4. Install Rnmr1D from GitHub using devtools
RUN R -e "devtools::install_github('INRA/Rnmr1D', dependencies = TRUE)"

# 5. Copy your Plumber API script into the container's filesystem
COPY plumber.R /app/plumber.R

# 6. Set the working directory inside the container
WORKDIR /app

# 7. Expose the port that the API will run on.
EXPOSE 8000

# 8. The command to run when the container starts. This starts the Plumber API server.
#    Render will automatically use the PORT environment variable if set.
CMD ["R", "-e", "port <- Sys.getenv('PORT', 8000); pr <- plumber::plumb('plumber.R'); pr$run(host='0.0.0.0', port=as.numeric(port))"]