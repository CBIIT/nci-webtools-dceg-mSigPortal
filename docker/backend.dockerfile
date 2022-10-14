FROM public.ecr.aws/amazonlinux/amazonlinux:2022

RUN dnf -y update \
    && dnf -y install \
    dnf-plugins-core \
    cairo-devel \
    cmake \
    git \
    gmp-devel \
    libcurl-devel \
    libjpeg-turbo-devel \
    libxml2-devel \
    mpfr-devel \
    nodejs \
    npm \
    openssl-devel \
    python3-devel \
    python3-pip \
    python3-setuptools \
    rsync \
    wget \ 
    which \
    R-4.1.3 \
    # && dnf -y builddep R \
    && dnf clean all

#  install R from source 
# ENV R_VERSION=4.1.3 
# RUN cd /tmp && \
#     curl -O https://cran.rstudio.com/src/base/R-4/R-${R_VERSION}.tar.gz && \
#     tar -xzvf R-${R_VERSION}.tar.gz 
# RUN cd /tmp/R-${R_VERSION} && \
#     ./configure \
#     --prefix=/opt/R/${R_VERSION} \
#     --enable-memory-profiling \
#     --enable-R-shlib \
#     --with-blas \
#     --with-lapack && \
#     make && \
#     make install
# RUN ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R && \
#     ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

RUN mkdir -p /deploy/server /deploy/logs

# install system fonts
RUN cd /tmp && \
    git clone https://github.com/xtmgah/SigProfilerPlotting && \
    cp /tmp/SigProfilerPlotting/fonts/* /usr/share/fonts && \
    fc-cache -fv;

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
# RUN pip3 install pandas 

# install client python packages
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerClusters#egg=SigProfilerClusters'
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerPlotting#egg=SigProfilerPlotting'
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'

# install genomes
# RUN python3.6 -c "\
# from SigProfilerMatrixGenerator import install as genInstall; \
# genInstall.install('GRCh37', rsync=False, bash=True); \
# genInstall.install('GRCh38', rsync=False, bash=True); \
# genInstall.install('mm10', rsync=False, bash=True)"



# use build cache for npm packages
COPY server/package*.json /deploy/server/

RUN npm install

# copy the rest of the application
COPY server .

CMD npm start

# docker build -t msigportal-backend -f backend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8330:8330 -v ~/Projects/msigportal/logs/:/deploy/logs -v ~/Projects/msigportal/tmp:/deploy/tmp -v ~/Projects/msigportal/config:/deploy/config -v ~/Projects/sigprofiler/data/genomes:/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb  -v ~/.aws/credentials:/root/.aws/credentials:ro --name msigportal-backend msigportal-backend 