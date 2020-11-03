ARG BASE_IMAGE=msigportal:base

FROM ${BASE_IMAGE}

# ARG SPARRPOWR_TAG

COPY . /deploy

# always install latest version of spatstat.core
# RUN Rscript -e "remotes::install_github('spatstat/spatstat.core')"

# install version of sparrpowR specified by tag
# RUN Rscript -e "remotes::install_github('machiela-lab/sparrpowR', ref='$SPARRPOWR_TAG')"

WORKDIR /deploy/server

RUN npm install

CMD npm start