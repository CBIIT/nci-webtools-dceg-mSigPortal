FROM public.ecr.aws/amazonlinux/amazonlinux:2023

RUN dnf -y update \
    && dnf -y install \
    bzip2 \
    cairo-devel \
    cmake \
    dnf-plugins-core \
    gcc \
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
    rsync \
    tar \
    wget \ 
    which \
    && dnf clean all

RUN mkdir -p /deploy/app /deploy/logs


# install Python packages from requirements.txt
WORKDIR /deploy/app
COPY extraction-service/requirements.txt /deploy/app/
RUN pip3 install -r requirements.txt


# install genomes
## NOTE: genomes do not need to be installed. They are saved on the host in [app]/data and mounted as a volume to the 
## sigprofilermatrixgenerator install directory. Verify path with "pip3 list"
# RUN python3 -c "from SigProfilerMatrixGenerator import install as genInstall; \
#     genInstall.install('GRCh37', rsync=False, bash=True)"
# genInstall.install('GRCh38', rsync=False, bash=True); \
# genInstall.install('mm10', rsync=False, bash=True)"


# use build cache for npm packages
COPY extraction-service/package*.json /deploy/app
RUN npm install

# copy the rest of the application
COPY extraction-service /deploy/app/

# ensure symlink exists for /data/genomes
ENV GENOME_PATH=/deploy/app/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb
RUN mkdir -p /data/genomes ${GENOME_PATH} \
    && rm -rf ${GENOME_PATH} \
    && ln -sf /data/genomes ${GENOME_PATH}

CMD npm run start
