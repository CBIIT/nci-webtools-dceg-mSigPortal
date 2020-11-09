ARG BASE_IMAGE=msigportal:base

FROM ${BASE_IMAGE}

COPY . /deploy

WORKDIR /deploy/server

RUN npm install

CMD npm start