FROM public.ecr.aws/amazonlinux/amazonlinux:2023

RUN dnf -y update \
    && dnf -y install \
    bzip2 \
    cairo-devel \
    cmake \
    dnf-plugins-core \
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
    python3-wheel \
    R-4.1.3 \
    rsync \
    tar \
    wget \ 
    which \
    && dnf clean all

# install nlopt
ENV LD_LIBRARY_PATH=/usr/local/lib64/
ENV NLOPT_VERSION=2.7.1
RUN cd /tmp \
    && curl -L https://github.com/stevengj/nlopt/archive/v${NLOPT_VERSION}.tar.gz | tar xz \
    && cd nlopt-${NLOPT_VERSION} \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install

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


# install python packages
RUN pip3 install pandas==1.5.3 seaborn 

# install client python packages
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerClusters#egg=SigProfilerClusters'
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerPlotting#egg=SigProfilerPlotting'
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'

# install bcftools
RUN cd /tmp \
    && curl -L https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 | tar xj \
    && cd bcftools-1.16 \
    && ./configure --enable-libcurl --prefix=/tmp/bcftools-1.16  \
    && make \
    && make install \
    && mv ./bcftools /usr/local/bin \
    && chmod +x /usr/local/bin/bcftools

# install genomes
## NOTE: genomes do not need to be installed. They are saved on the host in [app]/data and mounted as a volume to the 
## sigprofilermatrixgenerator install directory. Verify path with "pip3 list"
# RUN python3 -c "\
# from SigProfilerMatrixGenerator import install as genInstall; \
# genInstall.install('GRCh37', rsync=False, bash=True); \
# genInstall.install('GRCh38', rsync=False, bash=True); \
# genInstall.install('mm10', rsync=False, bash=True)"

# install R packages with renv
COPY server/renv.lock /deploy/server/
COPY server/.Rprofile /deploy/server/
COPY server/renv/activate.R /deploy/server/renv/
COPY server/renv/settings.dcf /deploy/server/renv/
COPY server/r-packages /deploy/server/r-packages

# copy renv cache if available
RUN mkdir /deploy/server/renv/.cache
ENV RENV_PATHS_CACHE=/deploy/server/renv/.cache
ARG R_RENV_CACHE_HOST=/renvCach[e]
COPY ${R_RENV_CACHE_HOST} ${RENV_PATHS_CACHE}
WORKDIR /deploy/server
RUN R -e "\
    options(\
    renv.config.repos.override = 'https://packagemanager.posit.co/cran/__linux__/rhel9/latest', \
    Ncpus = parallel::detectCores() \
    ); \
    renv::restore();"

# use build cache for npm packages
COPY server/package*.json /deploy/server/

RUN npm install

# copy the rest of the application
ARG CACHE_BUST
COPY server .

CMD npm start

# ensure symlink exists for /data/genomes
ENV GENOME_PATH=/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb
RUN mkdir -p /data/genomes ${GENOME_PATH} \
    && rm -rf ${GENOME_PATH} \
    && ln -sf /data/genomes ${GENOME_PATH}

# docker build -t msigportal-backend -f backend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8330:8330 -v ~/Projects/msigportal/logs/:/deploy/logs -v ~/Projects/msigportal/tmp:/deploy/tmp -v ~/Projects/msigportal/config:/deploy/config -v ~/Projects/sigprofiler/data/genomes:/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb  -v ~/.aws/credentials:/root/.aws/credentials:ro --name msigportal-backend msigportal-backend 