ARG BASE_IMAGE=msigportal:base

FROM ${BASE_IMAGE}

COPY . /deploy

RUN unzip /deploy/data/Database.zip -d /deploy/data/

WORKDIR /deploy/server

RUN npm install

CMD npm start