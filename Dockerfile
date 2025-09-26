# 1. Start from a trusted, official R base image
FROM rocker/r-ver:4.3.2

# 2. Install system dependencies required by R packages.
#    The XML package (a dependency of Rnmr1D) needs libxml2-dev.
RUN apt-get update -qq && apt-get install -y \
  libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# 3. Install the R packages you need from CRAN.
#    This step can take a few minutes during the first build.
RUN R -e "install.packages(c('plumber', 'Rnmr1D'))"

# 4. Copy your Plumber API script into the container's filesystem
COPY plumber.R /app/plumber.R

# 5. Set the working directory inside the container
WORKDIR /app

# 6. Expose the port that the API will run on. DigitalOcean needs this.
EXPOSE 8000

# 7. The command to run when the container starts. This starts the Plumber API server.
#    host='0.0.0.0' is required to make it accessible from outside the container.
CMD ["R", "-e", "pr <- plumber::plumb('plumber.R'); pr$run(host='0.0.0.0', port=8000)"]