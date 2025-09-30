# 1. Start from a trusted, official R base image
FROM rocker/r-ver:4.3.2

# 2. Install ALL required system-level dependencies.
#    - libnetcdf-dev, libhdf5-dev, zlib1g-dev: REQUIRED for mzR package to read mass-spec files.
#    - libglpk-dev: The GNU Linear Programming Kit, for 'igraph' (dependency of another package).
#    - libfftw3-dev: For Fast Fourier Transform, required by Rnmr1D.
#    - libarchive-dev: For the R 'archive' package.
#    - Others: Common build-time dependencies for R packages.
RUN apt-get update -qq && apt-get install -y \
    # For mzR
    libnetcdf-dev \
    libhdf5-dev \
    zlib1g-dev \
    # For other packages
    libarchive-dev \
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
        "ggplot2", "base64enc","plotly","dplyr"                          \
    ));                                                                  \
    BiocManager::install(c(                                              \
        "impute", "MassSpecWavelet", "pcaMethods",                       \
        "Spectra","mzR"                                                  \
    ), update=FALSE);                                                    \
'

# 4. Install Rnmr1D from GitHub using devtools
RUN R -e "devtools::install_github('INRA/Rnmr1D', dependencies = TRUE)"

# 5. Copy your Plumber API script into the container's filesystem
COPY plumber.R /app/plumber.R
COPY scripts/* /app/scripts/

# 6. Set the working directory inside the container
WORKDIR /app

# 7. Expose the port that the API will run on.
EXPOSE 8000

# 8. The command to run when the container starts. This starts the Plumber API server.
#    Render will automatically use the PORT environment variable if set.
CMD ["R", "-e", "port <- Sys.getenv('PORT', 8000); pr <- plumber::plumb('plumber.R'); pr$run(host='0.0.0.0', port=as.numeric(port))"]