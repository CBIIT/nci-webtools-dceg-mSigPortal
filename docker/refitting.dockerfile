FROM public.ecr.aws/amazonlinux/amazonlinux:2023

RUN dnf -y update \
    && dnf -y install \
    nodejs20 \
    nodejs20-npm  \
    tar \ 
    gzip \
    libcurl-devel \
    wget \
    R \
    && dnf clean all

RUN ln -s -f /usr/bin/node-20 /usr/bin/node; ln -s -f /usr/bin/npm-20 /usr/bin/npm;
RUN mkdir -p /refitting-service

WORKDIR /refitting-service

# Install essential R packages manually (minimal approach)
RUN R --vanilla -e "\
    # Configure repositories \
    options(\
        repos = c(CRAN = 'https://cloud.r-project.org'),\
        Ncpus = parallel::detectCores(),\
        timeout = 600\
    ); \
    # Install BiocManager first \
    install.packages('BiocManager', dependencies = FALSE);\
    library(BiocManager);\
    # Set compatible Bioconductor version \
    BiocManager::install(version = '3.18', ask = FALSE, update = FALSE);\
    # Install essential packages in stages \
    # Stage 1: Basic packages \
    install.packages(c('dplyr', 'ggplot2', 'readr', 'tibble'), dependencies = TRUE);\
    # Stage 2: Core Bioconductor packages \
    BiocManager::install(c('BiocGenerics', 'S4Vectors', 'IRanges'), ask = FALSE, update = FALSE);\
    # Stage 3: Advanced Bioconductor packages \
    BiocManager::install(c('GenomicRanges', 'Biostrings'), ask = FALSE, update = FALSE);\
    # Stage 4: BSgenome (most complex package) \
    tryCatch({\
        BiocManager::install('BSgenome', ask = FALSE, update = FALSE);\
    }, error = function(e) {\
        message('BSgenome installation failed, will continue without it: ', e\$message);\
    });\
    message('R package installation completed');\
    "

# Copy renv files for reference (but don't restore yet)
COPY refitting-service/renv.lock /refitting-service/
COPY refitting-service/.Rprofile /refitting-service/
COPY refitting-service/renv/activate.R /refitting-service/renv/
COPY refitting-service/renv/settings.json /refitting-service/renv/

COPY refitting-service /refitting-service/

COPY refitting-service/package.json refitting-service/package-lock.json ./

RUN npm install

# Create ENV file if it doesn't exist https://github.com/nodejs/node/issues/50993
RUN touch .env

# CMD npm start
CMD ["sleep", "infinity"]
