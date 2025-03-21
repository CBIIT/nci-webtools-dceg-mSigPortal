FROM nvidia/cuda:12.1.1-runtime-rockylinux9

RUN dnf -y update \
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


# Add CUDA repo and install cuSPARSE Lt to get libcusparseLt.so.0
RUN dnf config-manager --add-repo https://developer.download.nvidia.com/compute/cuda/repos/rhel9/x86_64/cuda-rhel9.repo \
    && dnf clean all \
    && dnf -y install cuda-cusparse-12-1 cuda-cusparse-devel-12-1


RUN mkdir -p /deploy/app /deploy/logs
WORKDIR /deploy/app


# install pytorch
RUN pip3 install torch
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/python3.9/site-packages/nvidia/cudnn/lib/:/usr/local/lib/python3.9/site-packages/nvidia/cuda_cupti/lib/

# install SigProfilerExtractor
RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerExtractor#egg=SigProfilerExtractor'

# install other python packages
#RUN pip3 install pandas==1.5.3 PyPDF2==2.11.2 SigProfilerAssignment==0.0.14 sigProfilerPlotting==1.2.2 
RUN pip3 install --no-cache-dir PyPDF2==2.11.2 SigProfilerAssignment==0.0.14 sigProfilerPlotting==1.2.2
RUN pip3 install -e 'git+https://github.com/AlexandrovLab/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'


# Force uninstall any existing NumPy/Pandas versions
RUN pip3 uninstall -y numpy pandas \
    && pip3 install --no-cache-dir numpy==1.26.4 pandas==1.3.5

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

# ensure symlink exists for /data/genomes
ENV GENOME_PATH=/deploy/app/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb
RUN mkdir -p /data/genomes ${GENOME_PATH} \
    && rm -rf ${GENOME_PATH} \
    && ln -sf /data/genomes ${GENOME_PATH}

CMD npm run start
