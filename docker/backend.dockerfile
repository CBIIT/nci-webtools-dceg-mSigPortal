FROM ${BACKEND_BASE_IMAGE:-oraclelinux:8-slim}

RUN rpm -i https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm && \
    microdnf -y update && \
    microdnf -y module enable nodejs:14 && \
    microdnf -y --enablerepo=ol8_codeready_builder install \
        # NLopt \
        # NLopt-devel \
        cairo \
        cairo-devel \
        git \
        gmp-devel \
        google-roboto-condensed-fonts \
        libcurl-devel \
        libjpeg-turbo-devel \
        libxml2-devel \
        mpfr-devel \
        nodejs \
        npm \
        openssl-devel \
        python3-devel \
        python3-pip \
        R \
        rsync \
        wget \
    && microdnf clean all

# configure C++ Toolchain for installing dependency RStan - https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
ENV MAKEFLAGS='-j2'
RUN mkdir -p $HOME/.R && \
    echo "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC","CXX14=g++" >> $HOME/.R/Makevars

# install R packages
RUN Rscript -e "install.packages(\
    c(\
        'aws.ec2metadata', \
        'aws.s3', \
        'BiocManager', \
        'car', \
        'coop', \
        'cowplot', \
        'entropy', \
        'factoextra', \
        'FactoMineR', \
        'ggdendro', \
        'ggExtra', \
        'ggforce', \
        'ggpubr', \
        'ggrepel', \
        'ggridges', \
        'ggsci', \
        'ggside', \
        'ggstatsplot', \
        'ggtext', \
        'hrbrthemes', \
        'janitor', \
        'scales', \
        'svglite', \
        'tidyverse' \
    ), \
    repos='https://cloud.r-project.org/', \
    Ncpus = 2)"

# install Bioc Packages
RUN Rscript -e "BiocManager::install(\
    c(\
        'BSgenome.Hsapiens.UCSC.hg19', \
        'BSgenome.Hsapiens.UCSC.hg38', \
        'ChIPseeker', \
        'org.Hs.eg.db', \
        'TxDb.Hsapiens.UCSC.hg19.knownGene', \
        'TxDb.Hsapiens.UCSC.hg38.knownGene' \
    ))"

# install python packages
RUN pip3 install scipy statsmodels

# install system fonts
RUN cd /tmp && \
    git clone https://github.com/xtmgah/SigProfilerPlotting && \
    cp /tmp/SigProfilerPlotting/fonts/* /usr/share/fonts && \
    fc-cache -fv;

# install client python packages
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'

RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerPlotting#egg=SigProfilerPlotting'

# install genomes
# RUN python3.6 -c "\
# from SigProfilerMatrixGenerator import install as genInstall; \
# genInstall.install('GRCh37', rsync=False, bash=True); \
# genInstall.install('GRCh38', rsync=False, bash=True); \
# genInstall.install('mm10', rsync=False, bash=True)"

RUN mkdir -p /deploy/server /deploy/logs

WORKDIR /deploy/server

# use build cache for npm packages
COPY server/package*.json /deploy/server/

RUN npm install

# copy the rest of the application
COPY . /deploy/

CMD npm start

# docker build -t msigportal-backend -f backend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8330:8330 -v ~/Projects/msigportal/logs/:/deploy/logs -v ~/Projects/msigportal/tmp:/deploy/tmp -v ~/Projects/msigportal/config:/deploy/config -v ~/Projects/sigprofiler/tsb:/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb -e AWS_ACCESS_KEY_ID=key -e AWS_SECRET_ACCESS_KEY=key --name msigportal-backend msigportal-backend 