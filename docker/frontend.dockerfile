FROM public.ecr.aws/amazonlinux/amazonlinux:2023

RUN dnf -y update \
    && dnf -y install \
    gcc-c++ \
    httpd \
    make \
    nodejs \
    npm \
    && dnf clean all

RUN mkdir -p /app/client

WORKDIR /app/client

COPY client/package.json /app/client/

RUN npm install

ARG CACHE_BUST
ARG RELEASE_VERSION
ARG LAST_UPDATED

COPY client /app/client/

# Set environment variables for the build
ENV VITE_APP_VERSION=${RELEASE_VERSION:-1.0.0}
ENV VITE_APP_LAST_UPDATE=${LAST_UPDATED:-"2025-09-04"}

RUN npm run build

RUN mkdir -p /var/www/html/mutational-signatures \
    && cp -r /app/client/dist/* /var/www/html/mutational-signatures

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


# docker build -t msigportal-frontend -f frontend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8331:80 --name msigportal-frontend msigportal-frontend