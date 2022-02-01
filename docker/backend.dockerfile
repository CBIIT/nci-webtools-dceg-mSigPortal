FROM ${BACKEND_BASE_IMAGE:-quay.io/centos/centos:stream8}

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
    cmake \
    v8-devel \
    # https://download.fedoraproject.org/pub/epel/7/x86_64/Packages/j/jq-1.6-2.el7.x86_64.rpm \
    # https://download.fedoraproject.org/pub/epel/7/x86_64/Packages/j/jq-devel-1.6-2.el7.x86_64.rpm \
    libjpeg-turbo-devel \
    openssl-devel \
    nodejs \
    R \
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
    # NLopt \
    # NLopt-devel \
    # google-roboto-condensed-fonts \
    && dnf clean all

# configure C++ Toolchain for installing dependency RStan - https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
ENV MAKEFLAGS='-j2'
RUN mkdir -p $HOME/.R && \
    echo -e "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC \nCXX14=g++" >> $HOME/.R/Makevars

RUN mkdir -p /deploy/server /deploy/logs

WORKDIR /deploy/server

COPY server/install.R .

RUN Rscript install.R

# install python packages
RUN pip3 install scipy statsmodels

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

# use build cache for npm packages
COPY server/package*.json /deploy/server/

RUN npm install

# copy the rest of the application
COPY . /deploy/

CMD npm start

# docker build -t msigportal-backend -f backend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8330:8330 -v ~/Projects/msigportal/logs/:/deploy/logs -v ~/Projects/msigportal/tmp:/deploy/tmp -v ~/Projects/msigportal/config:/deploy/config -v ~/Projects/sigprofiler/tsb:/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb  -v ~/.aws:/root/.aws:ro --name msigportal-backend msigportal-backend 