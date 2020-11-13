ARG BASE_IMAGE=msigportal:base

FROM ${BASE_IMAGE}

COPY server/ /deploy/server

WORKDIR /deploy/server

RUN npm install

CMD npm start