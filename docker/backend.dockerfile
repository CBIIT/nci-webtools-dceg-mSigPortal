FROM public.ecr.aws/amazonlinux/amazonlinux:2023

RUN dnf -y update \
    && dnf -y install spal-release \
    && dnf -y install \
    dnf-plugins-core \
    git \
    tar \
    nodejs24 \
    python3.13 \
    python3.13-devel \
    python3.13-pip \
    R-4.3.2 \
    bcftools \
    glpk \
    libuv \
    libcurl-devel \
    && dnf clean all

RUN mkdir -p /deploy/server /deploy/logs

# install system fonts
RUN cd /tmp && \
    git clone https://github.com/xtmgah/SigProfilerPlotting && \
    cp /tmp/SigProfilerPlotting/fonts/* /usr/share/fonts && \
    fc-cache -fv;

# install python packages (system python3 stays at 3.9; the app runs $PYTHON_BIN)
ENV PYTHON_BIN=python3.13
ARG PY_PINS="numpy>=2.1,<3 pandas>=2.2.3,<3 matplotlib>=3.9 seaborn==0.13.2"
RUN python3.13 -m pip install --no-cache-dir -e 'git+https://github.com/xtmgah/SigProfilerPlotting@python313#egg=SigProfilerPlotting' ${PY_PINS}
RUN python3.13 -m pip install --no-cache-dir -e 'git+https://github.com/xtmgah/SigProfilerMatrixGenerator@python313#egg=SigProfilerMatrixGenerator' ${PY_PINS}
# Simulator and Clusters declare PyPI SPMG/SPP (and Extractor) deps that would replace the forks
RUN python3.13 -m pip install --no-cache-dir --no-deps SigProfilerSimulator==1.2.2 \
    && python3.13 -m pip install --no-cache-dir 'fastrand>=1.2'
RUN python3.13 -m pip install --no-cache-dir --no-deps -e 'git+https://github.com/xtmgah/SigProfilerClusters@python313#egg=SigProfilerClusters' \
    && python3.13 -m pip install --no-cache-dir statsmodels scipy ${PY_PINS}

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

RUN mkdir -p /root/.R && echo "CXX11STD = -std=gnu++14" >> /root/.R/Makevars
RUN Rscript install.R

# install npm packages
COPY server/package*.json /deploy/server/

RUN npm install

# copy the rest of the application
COPY server .

# Create ENV file so nodejs doesn't crash https://github.com/nodejs/node/issues/50993
RUN touch .env

# ensure symlink exists for /data/genomes
ENV GENOME_PATH=/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb
RUN mkdir -p /data/genomes ${GENOME_PATH} \
    && rm -rf ${GENOME_PATH} \
    && ln -sf /data/genomes ${GENOME_PATH}

CMD npm start