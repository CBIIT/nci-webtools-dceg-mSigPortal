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
    R-4.3.2 \
    fribidi-devel \
    libtiff-devel \
    rsync \
    tar \
    wget \ 
    which \
    && dnf clean all


RUN mkdir -p /deploy/server /deploy/logs

# install system fonts
RUN cd /tmp && \
    git clone https://github.com/xtmgah/SigProfilerPlotting && \
    cp /tmp/SigProfilerPlotting/fonts/* /usr/share/fonts && \
    fc-cache -fv;


# install python packages
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerClusters#egg=SigProfilerClusters'
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerPlotting#egg=SigProfilerPlotting'
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'
RUN pip3 install --force-reinstall --no-cache-dir numpy==1.26.4 pandas==1.3.5

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
# NOTE: genomes do not need to be installed. They are saved on the host in [app]/data and mounted as a volume to the 
# sigprofilermatrixgenerator install directory. Verify path with "pip3 list"
# RUN python3 -c "\
# from SigProfilerMatrixGenerator import install as genInstall; \
# genInstall.install('GRCh37', rsync=False, bash=True); \
# genInstall.install('GRCh38', rsync=False, bash=True); \
# genInstall.install('mm10', rsync=False, bash=True)"

WORKDIR /deploy/server

# install R packages
COPY server/.Rprofile .
COPY server/install.R .
COPY server/r-packages .

RUN Rscript install.R

# install npm packages
COPY server/package*.json /deploy/server/

RUN npm install

# copy the rest of the application
COPY server .


# ensure symlink exists for /data/genomes
ENV GENOME_PATH=/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb
RUN mkdir -p /data/genomes ${GENOME_PATH} \
    && rm -rf ${GENOME_PATH} \
    && ln -sf /data/genomes ${GENOME_PATH}

CMD npm start