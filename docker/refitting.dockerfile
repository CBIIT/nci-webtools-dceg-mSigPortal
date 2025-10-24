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
# Install R packages in separate steps for better error handling
# Step 1: Install BiocManager and basic setup
RUN R --vanilla -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org'), timeout = 600); \
    if (!('BiocManager' %in% rownames(installed.packages()))) { \
        install.packages('BiocManager', dependencies = FALSE); \
    }; \
    library(BiocManager); \
    BiocManager::install(version = '3.21', ask = FALSE, update = FALSE); \
    message('BiocManager setup completed'); \
    "

# Step 2: Install tidyverse (most common dependency issue)
RUN R --vanilla -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org'), timeout = 600); \
    if (!('tidyverse' %in% rownames(installed.packages()))) { \
        message('Installing tidyverse...'); \
        install.packages('tidyverse', dependencies = TRUE); \
    } else { \
        message('tidyverse already installed'); \
    }; \
    message('Tidyverse installation completed'); \
    "

# Step 3: Install basic CRAN packages including jsonlite
RUN R --vanilla -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org'), timeout = 600); \
    packages <- c('dplyr', 'ggplot2', 'readr', 'tibble', 'jsonlite'); \
    for (pkg in packages) { \
        if (!(pkg %in% rownames(installed.packages()))) { \
            message('Installing: ', pkg); \
            install.packages(pkg, dependencies = TRUE); \
        } else { \
            message('Already installed: ', pkg); \
        } \
    }; \
    message('Basic CRAN packages completed'); \
    "

# Step 4: Install Bioconductor core packages
RUN R --vanilla -e "\
    library(BiocManager); \
    packages <- c('BiocGenerics', 'S4Vectors', 'IRanges'); \
    for (pkg in packages) { \
        if (!(pkg %in% rownames(installed.packages()))) { \
            message('Installing Bioc package: ', pkg); \
            BiocManager::install(pkg, ask = FALSE, update = FALSE); \
        } else { \
            message('Already installed: ', pkg); \
        } \
    }; \
    message('Bioconductor core packages completed'); \
    "

# Step 5: Install advanced Bioconductor packages
RUN R --vanilla -e "\
    library(BiocManager); \
    packages <- c('GenomicRanges', 'Biostrings'); \
    for (pkg in packages) { \
        if (!(pkg %in% rownames(installed.packages()))) { \
            message('Installing advanced Bioc package: ', pkg); \
            BiocManager::install(pkg, ask = FALSE, update = FALSE); \
        } else { \
            message('Already installed: ', pkg); \
        } \
    }; \
    message('Advanced Bioconductor packages completed'); \
    "

# Step 6: Install SATS (with error handling)
RUN R --vanilla -e "\
    tryCatch({ \
        if (!('SATS' %in% rownames(installed.packages()))) { \
            message('Installing SATS...'); \
            install.packages('SATS', dependencies = TRUE); \
        } else { \
            message('SATS already installed'); \
        } \
    }, error = function(e) { \
        message('SATS installation failed, will continue: ', e\$message); \
    }); \
    message('SATS installation step completed'); \
    "

# Step 7: Install BSgenome (with error handling)
RUN R --vanilla -e "\
    library(BiocManager); \
    tryCatch({ \
        if (!('BSgenome' %in% rownames(installed.packages()))) { \
            message('Installing BSgenome...'); \
            BiocManager::install('BSgenome', ask = FALSE, update = FALSE); \
        } else { \
            message('BSgenome already installed'); \
        } \
    }, error = function(e) { \
        message('BSgenome installation failed, will continue: ', e\$message); \
    }); \
    message('BSgenome installation step completed'); \
    "

# Step 8: Install BSgenome reference packages required by SATS
RUN R --vanilla -e "\
    library(BiocManager); \
    packages <- c('BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'); \
    for (pkg in packages) { \
        tryCatch({ \
            if (!(pkg %in% rownames(installed.packages()))) { \
                message('Installing BSgenome reference: ', pkg); \
                BiocManager::install(pkg, ask = FALSE, update = FALSE); \
            } else { \
                message('Already installed: ', pkg); \
            } \
        }, error = function(e) { \
            message('Failed to install ', pkg, ': ', e\$message); \
        }); \
    }; \
    message('BSgenome reference packages completed'); \
    "

# Step 9: Verify all essential packages are available
RUN R --vanilla -e "\
    essential_packages <- c('tidyverse', 'GenomicRanges', 'Biostrings', 'SATS', 'BSgenome', 'jsonlite'); \
    missing_packages <- c(); \
    for (pkg in essential_packages) { \
        if (!(pkg %in% rownames(installed.packages()))) { \
            missing_packages <- c(missing_packages, pkg); \
        } \
    }; \
    if (length(missing_packages) > 0) { \
        message('WARNING: Missing essential packages: ', paste(missing_packages, collapse=', ')); \
    } else { \
        message('SUCCESS: All essential packages are installed'); \
    }; \
    "

# Copy renv files for reference (but don't restore yet)
COPY refitting-service/renv.lock /refitting-service/
COPY refitting-service/.Rprofile /refitting-service/
COPY refitting-service/renv/activate.R /refitting-service/renv/
COPY refitting-service/renv/settings.json /refitting-service/renv/

COPY refitting-service /refitting-service/

# Copy package files and install npm dependencies
COPY refitting-service/package.json refitting-service/package-lock.json ./

RUN npm install

# Create ENV file if it doesn't exist https://github.com/nodejs/node/issues/50993
RUN touch .env

# Test that the refitting R script loads correctly
RUN R --vanilla -e "\
    tryCatch({ \
        source('./refitting.R'); \
        if (exists('run_sbs_refitting')) { \
            message('SUCCESS: refitting.R loaded and run_sbs_refitting function is available'); \
        } else { \
            stop('ERROR: run_sbs_refitting function not found'); \
        } \
    }, error = function(e) { \
        message('ERROR loading refitting.R: ', e\$message); \
        quit(status=1); \
    }); \
    "

# Expose the port that the refitting service uses
EXPOSE 8334

# Start the refitting service
CMD ["npm", "start"]
