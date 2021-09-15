FROM quay.io/centos/centos:stream8

RUN dnf -y update \
    && dnf -y install \
    dnf-plugins-core \
    epel-release \
    glibc-langpack-en \
    && dnf config-manager --enable powertools \
    && dnf -y module enable nodejs:13 \
    && dnf -y install \
    # gdal-devel \
    # proj-devel \
    # protobuf-devel \
    # udunits2-devel \
    v8-devel \
    # https://download.fedoraproject.org/pub/epel/7/x86_64/Packages/j/jq-1.6-2.el7.x86_64.rpm \
    # https://download.fedoraproject.org/pub/epel/7/x86_64/Packages/j/jq-devel-1.6-2.el7.x86_64.rpm \
    libjpeg-turbo-devel \
    openssl-devel \
    nodejs \
    R \
    python3-pip \
    libcurl-devel \
    libxml2-devel \
    git \
    rsync \
    wget \
    && dnf -y install \
    gmp-devel \
    mpfr-devel \
    cairo \
    cairo-devel \
    # NLopt \
    # NLopt-devel \
    google-roboto-condensed-fonts \
    && dnf clean all

# install python packages
RUN pip3 install scipy statsmodels

# install system fonts
RUN cd /tmp; git clone https://github.com/xtmgah/SigProfilerPlotting; cp /tmp/SigProfilerPlotting/fonts/* /usr/share/fonts; fc-cache -fv;

# install client python packages
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'

RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerPlotting#egg=SigProfilerPlotting'

# install genomes
RUN python3.6 -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh37', rsync=False, bash=True); genInstall.install('GRCh38', rsync=False, bash=True); genInstall.install('mm10', rsync=False, bash=True)"

# configure C++ Toolchain for installing dependency RStan - https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
RUN Rscript -e 'dotR <- file.path(Sys.getenv("HOME"), ".R"); dir.create(dotR); M <- file.path(dotR, "Makevars"); file.create(M); cat("\nCXX14FLAGS=-O3 -march=native -mtune=native -fPIC","CXX14=g++", file = M, sep = "\n", append = TRUE);'

# install R packages
RUN Rscript -e "Sys.setenv(MAKEFLAGS = '-j2'); install.packages(c('tidyverse', 'hrbrthemes', 'ggExtra', 'ggsci', 'ggrepel', 'ggdendro', 'scales', 'ggforce', 'svglite', 'cowplot', 'car', 'FactoMineR', 'factoextra', 'coop', 'ggridges', 'ggtext', 'ggpubr', 'entropy', 'janitor', 'ggstatsplot', 'BiocManager', 'aws.s3', 'aws.ec2metadata'), repos='https://cloud.r-project.org/', Ncpus = 2)"

# install Bioc Packages
RUN Rscript -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene", "BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "ChIPseeker", "org.Hs.eg.db"))'

RUN mkdir -p /deploy/server /deploy/logs

WORKDIR /deploy/server

# use build cache for npm packages
COPY server/package*.json /deploy/server/

RUN npm install

# copy the rest of the application
COPY . /deploy/

CMD npm start