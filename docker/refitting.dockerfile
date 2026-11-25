FROM public.ecr.aws/amazonlinux/amazonlinux:2023

RUN dnf -y update \
    && dnf -y install \
    nodejs24 \
    R-4.3.2 \
    tar \ 
    gzip \
    libcurl-devel \
    && dnf clean all

RUN mkdir -p /refitting-service

WORKDIR /refitting-service

COPY refitting-service/.Rprofile .
COPY refitting-service/install.R .

RUN Rscript install.R

COPY refitting-service/package.json refitting-service/package-lock.json ./

RUN npm install

COPY refitting-service .

# Create ENV file if it doesn't exist https://github.com/nodejs/node/issues/50993
RUN touch .env

# Remove Python 3.9 and its dependents (including dnf/yum) after all package installs.
RUN dnf -y --setopt=protected_packages= --setopt=clean_requirements_on_remove=False remove python3

EXPOSE 8334

CMD npm run start
