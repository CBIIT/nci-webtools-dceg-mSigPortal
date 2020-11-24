FROM centos:latest

RUN dnf -y update \
   && dnf -y install \
   dnf-plugins-core \
   epel-release \
   glibc-langpack-en \
   && dnf config-manager --enable PowerTools \
   && dnf -y module enable nodejs:13 \
   && dnf -y install \
   # gdal-devel \
   # proj-devel \
   # protobuf-devel \
   # udunits2-devel \
   v8-devel \
   # https://download.fedoraproject.org/pub/epel/7/x86_64/Packages/j/jq-1.6-2.el7.x86_64.rpm \
   # https://download.fedoraproject.org/pub/epel/7/x86_64/Packages/j/jq-devel-1.6-2.el7.x86_64.rpm \
   libjpeg-turbo-devel \
   openssl-devel \
   nodejs \
   R \
   python3-pip \
   libcurl-devel \
   libxml2-devel \
   git \
   rsync \
   wget \
   && dnf -y install \
   gmp-devel \
   mpfr-devel \
   cairo \
   cairo-devel \
   # NLopt \
   # NLopt-devel \
   google-roboto-condensed-fonts \
   && dnf clean all

# COPY . /deploy/data

# RUN unzip /deploy/data/Database.zip -d /deploy/data/

RUN pip3 install scipy statsmodels

RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerMatrixGenerator#egg=SigProfilerMatrixGenerator'

RUN pip3 install -e 'git+https://github.com/xtmgah/SigProfilerPlotting#egg=SigProfilerPlotting'

RUN python3.6 -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh37', rsync=False, bash=True); genInstall.install('GRCh38', rsync=False, bash=True); genInstall.install('mm10', rsync=False, bash=True)"

RUN R -e "Sys.setenv(MAKEFLAGS = '-j2'); install.packages('ggstatsplot', repos='https://cloud.r-project.org/')"

RUN R -e "Sys.setenv(MAKEFLAGS = '-j2'); install.packages(c('tidyverse', 'hrbrthemes', 'ggsci', 'ggrepel', 'ggdendro', 'scales', 'ggforce', 'svglite', 'cowplot', 'car', 'FactoMineR', 'factoextra', 'coop', 'ggridges', 'ggtext', 'ggpubr', 'entropy', 'janitor'), repos='https://cloud.r-project.org/')"