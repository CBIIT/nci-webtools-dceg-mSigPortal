FROM public.ecr.aws/amazonlinux/amazonlinux:2022

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
    rsync \
    tar \
    wget \ 
    which \
    gcc \
    && dnf clean all

RUN mkdir -p /deploy/app /deploy/logs

# install SigProfilerExtractor from local zip
COPY packages/SigProfilerExtractor/ /deploy/app/SigProfilerExtractor
WORKDIR /deploy/app
RUN pip3 install -e SigProfilerExtractor/SigProfilerExtractor 

# install other python packages
RUN pip3 install SigProfilerAssignment==0.0.14 sigProfilerPlotting==1.2.2 PyPDF2==2.11.2
RUN pip3 install -e 'git+https://github.com/AlexandrovLab/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'


# install bcftools
# RUN cd /tmp \
#     && curl -L https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 | tar xj \
#     && cd bcftools-1.16 \
#     && ./configure --enable-libcurl --prefix=/tmp/bcftools-1.16  \
#     && make \
#     && make install \
#     && mv ./bcftools /usr/local/bin \
#     && chmod +x /usr/local/bin/bcftools

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
COPY extraction-service /deploy/app
COPY server /deploy

CMD npm run start-extraction-worker

# docker build -t msigportal-backend -f backend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8330:8330 -v ~/Projects/msigportal/logs/:/deploy/logs -v ~/Projects/msigportal/tmp:/deploy/tmp -v ~/Projects/msigportal/config:/deploy/config -v ~/Projects/sigprofiler/data/genomes:/src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb  -v ~/.aws/credentials:/root/.aws/credentials:ro --name msigportal-backend msigportal-backend 