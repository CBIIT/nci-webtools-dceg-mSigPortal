ARG BASE_IMAGE=msigportal:base

FROM ${BASE_IMAGE}

COPY . /deploy

RUN unzip /deploy/data/Database.zip

WORKDIR /deploy/server

RUN npm install

CMD npm start