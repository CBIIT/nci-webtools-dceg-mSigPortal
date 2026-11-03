FROM public.ecr.aws/amazonlinux/amazonlinux:2023 AS builder

RUN dnf -y update \
    && dnf -y install nodejs24 \
    && dnf clean all

WORKDIR /app/client

COPY client/package.json client/package-lock.json ./

RUN npm ci

COPY client/ ./

ARG VITE_APP_VERSION=local
ARG VITE_LAST_UPDATED
RUN echo "VITE_APP_VERSION=${VITE_APP_VERSION}" >> .env \
    && echo "VITE_LAST_UPDATED=${VITE_LAST_UPDATED}" >> .env

RUN npm run build

FROM public.ecr.aws/amazonlinux/amazonlinux:2023 AS runner

RUN dnf -y update \
    && dnf -y install httpd \
    && dnf clean all

COPY --from=builder /app/client/dist/ /var/www/html/mutational-signatures/

# Add custom httpd configuration
COPY docker/httpd-msigportal.conf /etc/httpd/conf.d/httpd-msigportal.conf

# forward request and error logs to docker log collector
RUN ln -sf /dev/stdout /var/log/httpd/access_log \
    && ln -sf /dev/stderr /var/log/httpd/error_log

WORKDIR /var/www/html

EXPOSE 80
EXPOSE 443

ENV SERVER_TIMEOUT=900

CMD rm -rf /run/httpd/* /tmp/httpd* \
    && exec /usr/sbin/httpd -DFOREGROUND
