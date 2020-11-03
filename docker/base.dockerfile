FROM centos:latest

RUN dnf -y update \
   && dnf -y install \
      dnf-plugins-core \
      epel-release \
      glibc-langpack-en \
   && dnf config-manager --enable PowerTools \
   && dnf -y module enable nodejs:13 \
   && dnf -y install \
      nodejs \
      R \
      python3-pip \
   && dnf -y install \
      cairo \
      cairo-devel \
      # NLopt \
      # NLopt-devel \
      google-roboto-condensed-fonts \
   && dnf clean all

RUN pip3 install scipy statsmodels

RUN Rscript -e "install.packages(c('hms', 'tidyverse', 'hrbrthemes', 'ggsci', 'ggrepel', 'ggdendro', 'ggscales', 'ggforce', 'svglite', 'cowplot', 'carData', 'pbkrtest', 'quantreg', 'car', 'FactoMineR', 'factoextra', 'coop', 'ggridges', 'ggstatsplot', 'ggtext'), repos='https://cloud.r-project.org/')"