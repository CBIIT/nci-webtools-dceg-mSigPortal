FROM alpine:latest

RUN apk upgrade && \
    apk add \
      automake \
      bash \
      cmake \
      freetype-dev \
      g++ \
      gfortran \
      git \
      lapack-dev \
      libffi-dev \
      libxml2-dev \
      libxslt-dev \
      nodejs \
      npm \
      py3-matplotlib \
      py3-pandas \
      py3-pip \
      py3-scipy \
      python3 \
      python3-dev \
      R \
      R-dev \
      subversion

# configure C++ Toolchain for installing dependency RStan - https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
RUN Rscript -e 'dotR <- file.path(Sys.getenv("HOME"), ".R"); dir.create(dotR); M <- file.path(dotR, "Makevars"); file.create(M); cat("\nCXX14FLAGS=-O3 -march=native -mtune=native -fPIC","CXX14=g++", file = M, sep = "\n", append = TRUE);'

# install R packages
RUN Rscript -e "Sys.setenv(MAKEFLAGS = '-j2'); install.packages(c('tidyverse', 'hrbrthemes', 'ggExtra', 'ggsci', 'ggrepel', 'ggdendro', 'scales', 'ggforce', 'svglite', 'cowplot', 'car', 'FactoMineR', 'factoextra', 'coop', 'ggridges', 'ggtext', 'ggpubr', 'entropy', 'janitor', 'ggstatsplot', 'BiocManager', 'aws.s3', 'aws.ec2metadata', 'ggside'), repos='https://cloud.r-project.org/', Ncpus = 2)"

# install Bioc Packages
RUN Rscript -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene", "BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "ChIPseeker", "org.Hs.eg.db"))'

# install python packages
RUN pip3 install statsmodels

# install system fonts from client python packages
RUN cd /tmp && \
    mkdir -p /usr/share/fonts && \
    git clone https://github.com/xtmgah/SigProfilerPlotting && \
    cp /tmp/SigProfilerPlotting/fonts/* /usr/share/fonts && \
    fc-cache -fv

# install client python packages
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerPlotting#egg=SigProfilerPlotting'

RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'

# install genomes
RUN python3 -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh37', rsync=False, bash=True); genInstall.install('GRCh38', rsync=False, bash=True); genInstall.install('mm10', rsync=False, bash=True)"

RUN mkdir -p /deploy/server /deploy/logs

WORKDIR /deploy/server

# use build cache for npm packages
COPY server/package*.json /deploy/server/

RUN npm install

# copy the rest of the application
COPY . /deploy/

CMD npm start