FROM nvidia/cuda:11.8.0-runtime-rockylinux8

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
    # python3-wheel \
    rsync \
    tar \
    wget \ 
    which \
    # https://us.download.nvidia.com/tesla/470.182.03/nvidia-driver-local-repo-rhel8-470.182.03-1.0-1.x86_64.rpm \
    && dnf clean all

RUN mkdir -p /deploy/app /deploy/logs

# install SigProfilerExtractor from local zip
COPY packages/SigProfilerExtractor/ /deploy/app/SigProfilerExtractor
WORKDIR /deploy/app
RUN pip3 install -e SigProfilerExtractor/SigProfilerExtractor 

# install other python packages
RUN pip3 install pandas==1.5.3 PyPDF2==2.11.2 SigProfilerAssignment==0.0.14 sigProfilerPlotting==1.2.2 
RUN pip3 install -e 'git+https://github.com/AlexandrovLab/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'


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
