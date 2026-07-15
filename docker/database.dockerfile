FROM public.ecr.aws/amazonlinux/amazonlinux:2023

RUN dnf -y update \
    && dnf -y install \
    dnf-plugins-core \
    nodejs24 \
    && dnf clean all

WORKDIR /deploy/database

COPY database/package*.json ./
RUN npm install --omit=dev

COPY database .

CMD ["node", "startDatabaseImport.js", "--schema", "schema.js", "--sources", "sources.js", "--provider", "s3"]
