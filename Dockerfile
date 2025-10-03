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
    curl \
    && rm -rf /var/lib/apt/lists/*

# 3. Install all required R packages from CRAN and Bioconductor in one go.
#    This is more efficient for Docker builds.
RUN R -e '                                                               \
    install.packages(c(                                                  \
        "archive",                                                       \
        "BiocManager", "devtools", "drc",                                \
        "ggplot2", "base64enc","dplyr"                                   \
    ));                                                                  \
    BiocManager::install(c(                                              \
        "impute", "MassSpecWavelet", "pcaMethods",                       \
        "Spectra","mzR"                                                  \
    ), update=FALSE);                                                    \
'

# 4. Install Rnmr1D from GitHub using devtools
RUN R -e "devtools::install_github('INRA/Rnmr1D', dependencies = TRUE)"
#    Install Node.js and npm using the official NodeSource repository
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
    && apt-get install -y nodejs
# 5. Copy your API script into the container's filesystem
WORKDIR /app
COPY server.js /app/server.js
COPY scripts/* /app/scripts/
COPY package*.json /app
RUN npm install
# 6. Set the working directory inside the container

# 7. Expose the port that the API will run on.
EXPOSE 8000

# 8. The command to run when the container starts. This starts the express API server.
CMD ["node", "server.js"]
