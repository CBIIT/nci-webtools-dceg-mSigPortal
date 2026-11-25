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

# Remove Python 3.9 and its dependents (including dnf/yum) after all package installs.
RUN dnf -y --setopt=protected_packages= --setopt=clean_requirements_on_remove=False remove python3

CMD ["node", "startDatabaseImport.js", "--schema", "schema.js", "--sources", "sources.js", "--provider", "s3"]
