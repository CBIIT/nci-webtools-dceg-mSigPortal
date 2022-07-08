FROM ${BACKEND_BASE_IMAGE:-quay.io/centos/centos:stream8}

ENV R_VERSION=4.1.2

RUN dnf -y update \
    && dnf -y install \
    dnf-plugins-core \
    epel-release \
    glibc-langpack-en \
    && dnf config-manager --enable powertools \
    && dnf -y module enable nodejs:13 \
    && dnf -y install \
    cmake \
    v8-devel \
    libjpeg-turbo-devel \
    openssl-devel \
    nodejs \
    python3-pip \
    python3-devel \
    libcurl-devel \
    libxml2-devel \
    git \
    rsync \
    wget \
    && dnf -y install \
    gmp-devel \
    mpfr-devel \
    cairo-devel \
    google-roboto-condensed-fonts \
    && dnf -y install https://cdn.rstudio.com/r/centos-8/pkgs/R-${R_VERSION}-1-1.x86_64.rpm \
    && dnf clean all

# symlink R
RUN ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
RUN ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

RUN mkdir -p /deploy/server /deploy/logs

# install system fonts
RUN cd /tmp && \
    git clone https://github.com/xtmgah/SigProfilerPlotting && \
    cp /tmp/SigProfilerPlotting/fonts/* /usr/share/fonts && \
    fc-cache -fv;

# install client python packages
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerPlotting#egg=SigProfilerPlotting'

RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'

# install genomes
# RUN python3.6 -c "\
# from SigProfilerMatrixGenerator import install as genInstall; \
# genInstall.install('GRCh37', rsync=False, bash=True); \
# genInstall.install('GRCh38', rsync=False, bash=True); \
# genInstall.install('mm10', rsync=False, bash=True)"

# configure C++ Toolchain for installing dependency RStan - https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
# ENV MAKEFLAGS='-j2'
# RUN mkdir -p $HOME/.R && \
#     echo -e "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC \nCXX14=g++" >> $HOME/.R/Makevars

# install renv
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org/')"

# install R packages
COPY server/renv.lock /deploy/server/

WORKDIR /deploy/server

RUN R -e "options(Ncpus=parallel::detectCores()); renv::restore()"

# install python packages
RUN pip3 install scipy statsmodels

# use build cache for npm packages
COPY server/package*.json /deploy/server/

RUN npm install

# copy the rest of the application
COPY server .

CMD npm start

# docker build -t msigportal-backend -f backend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8330:8330 -v ~/Projects/msigportal/logs/:/deploy/logs -v ~/Projects/msigportal/tmp:/deploy/tmp -v ~/Projects/msigportal/config:/deploy/config -v ~/Projects/sigprofiler/data/genomes:/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb  -v ~/.aws/credentials:/root/.aws/credentials:ro --name msigportal-backend msigportal-backend 