FROM public.ecr.aws/amazonlinux/amazonlinux:2023

RUN dnf -y update \
    && dnf -y install \
    nodejs20 \
    nodejs20-npm  \
    tar \ 
    gzip \
    libcurl-devel \
    wget \
    R-4.3.2 \
    && dnf clean all

RUN ln -s -f /usr/bin/node-20 /usr/bin/node; ln -s -f /usr/bin/npm-20 /usr/bin/npm;
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

EXPOSE 8334

CMD npm run start
