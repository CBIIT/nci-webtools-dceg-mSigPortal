FROM public.ecr.aws/amazonlinux/amazonlinux:2022

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
COPY client /app/client/

RUN npm run build

RUN mkdir -p /var/www/html/mutational-signatures \
    && cp -r /app/client/build/* /var/www/html/mutational-signatures

# Add custom httpd configuration
COPY docker/httpd-msigportal.conf /etc/httpd/conf.d/httpd-msigportal.conf

WORKDIR /var/www/html

EXPOSE 80
EXPOSE 443

CMD rm -rf /run/httpd/* /tmp/httpd* \
    && exec /usr/sbin/httpd -DFOREGROUND

# docker build -t msigportal-frontend -f frontend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8331:80 --name msigportal-frontend msigportal-frontend