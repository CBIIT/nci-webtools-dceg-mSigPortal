FROM public.ecr.aws/amazonlinux/amazonlinux:2023

RUN dnf -y update \
    && dnf -y install \
    nodejs20 \
    nodejs20-npm  \
    tar \ 
    R-4.3.2 \
    gzip \
    libcurl-devel \
    && dnf clean all

RUN ln -s -f /usr/bin/node-20 /usr/bin/node; ln -s -f /usr/bin/npm-20 /usr/bin/npm;
RUN mkdir -p /refitting-service

WORKDIR /refitting-service

# install R packages with renv
COPY refitting-service /refitting-service/
COPY refitting-service/.Rprofile /refitting-service/
COPY refitting-service/renv/activate.R /refitting-service/renv/
COPY refitting-service/renv/settings.json /refitting-service/renv/

RUN R -e "\
    options(\
    renv.config.repos.override = 'https://packagemanager.posit.co/cran/__linux__/rhel9/latest', \
    Ncpus = parallel::detectCores() \
    ); \
    renv::restore();"

COPY refitting-service/package.json refitting-service/package-lock.json ./

RUN npm install

# Create ENV file if it doesn't exist https://github.com/nodejs/node/issues/50993
RUN touch .env

# CMD npm start
CMD sleep infinity
