FROM nvidia/cuda:12.9.2-runtime-amzn2023

# the base is locked to an older AL2023 release that predates the python3.13 and nodejs24 packages
RUN dnf -y --releasever=latest update \
    && dnf -y install \
    bzip2 \
    cairo-devel \
    cmake \
    dnf-plugins-core \
    gcc \
    gcc-c++ \
    git \
    gmp-devel \
    libcurl-devel \
    libjpeg-turbo-devel \
    libxml2-devel \
    mpfr-devel \
    nodejs24 \
    openssl-devel \
    python3.13 \
    python3.13-devel \
    python3.13-pip \
    rsync \
    tar \
    wget \ 
    which \
    && dnf clean all


# Add CUDA repo and install cuSPARSE Lt to get libcusparseLt.so.0
# RUN dnf config-manager --add-repo https://developer.download.nvidia.com/compute/cuda/repos/rhel9/x86_64/cuda-rhel9.repo \
#     && dnf clean all \
#     && dnf -y install cuda-cusparse-12-1 cuda-cusparse-devel-12-1


RUN mkdir -p /deploy/app /deploy/logs
WORKDIR /deploy/app


# install pytorch (cu129 wheels to match the base image)
ENV PYTHON_BIN=python3.13
RUN python3.13 -m pip install --no-cache-dir "torch>=2.5,<3" --index-url https://download.pytorch.org/whl/cu129
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/python3.13/site-packages/nvidia/cudnn/lib/:/usr/local/lib/python3.13/site-packages/nvidia/cuda_cupti/lib/

# install SigProfilerExtractor (editable: a regular install would pick up the fork's stale committed build/lib)
RUN python3.13 -m pip install --no-cache-dir -e 'git+https://github.com/xtmgah/SigProfilerExtractor@python313#egg=SigProfilerExtractor' \
    "SigProfilerAssignment>=1.1.5,<1.2" "sigProfilerPlotting>=1.4.3,<1.5" "SigProfilerMatrixGenerator>=1.4,<1.5" \
    "pypdf>=6,<7" "reportlab>=4,<6" "numpy>=2.1,<3" "pandas>=2.2.3,<3" "matplotlib>=3.9,<4"

# install genomes
## NOTE: genomes do not need to be installed. They are saved on the host in [app]/data and mounted as a volume to the 
## sigprofilermatrixgenerator install directory. Verify path with "pip3 list"
# RUN python3 -c "from SigProfilerMatrixGenerator import install as genInstall; \
#     genInstall.install('GRCh37', rsync=False, bash=True)"
# genInstall.install('GRCh38', rsync=False, bash=True); \
# genInstall.install('mm10', rsync=False, bash=True)"


# use build cache for npm packages
COPY extraction-service/package*.json /deploy/app/
RUN npm install

# copy the rest of the application
COPY extraction-service /deploy/app/

# Create ENV file so nodejs doesn't crash https://github.com/nodejs/node/issues/50993
RUN touch .env

# ensure symlink exists for /data/genomes
ENV GENOME_PATH=/usr/local/lib/python3.13/site-packages/SigProfilerMatrixGenerator/references/chromosomes/tsb
RUN mkdir -p /data/genomes ${GENOME_PATH} \
    && rm -rf ${GENOME_PATH} \
    && ln -sf /data/genomes ${GENOME_PATH}

CMD npm run start
