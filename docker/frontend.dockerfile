FROM ${FRONTEND_BASE_IMAGE:-oraclelinux:8-slim}

RUN microdnf -y update \
    && microdnf -y module enable nodejs:14 \
    && microdnf -y install \
    gcc-c++ \
    httpd \
    make \
    nodejs \
    npm \
    && microdnf clean all

RUN mkdir /client

WORKDIR /client

COPY client/package*.json /client/

RUN npm install

COPY client /client/

RUN npm run build \
    && mv /client/build /var/www/html/mutational-signatures

# Add custom httpd configuration
COPY docker/httpd-msigportal.conf /etc/httpd/conf.d/httpd-msigportal.conf

WORKDIR /var/www/html

EXPOSE 80
EXPOSE 443

CMD rm -rf /run/httpd/* /tmp/httpd* \
    && exec /usr/sbin/apachectl -DFOREGROUND

# docker build -t msigportal-frontend -f frontend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8331:80 --name msigportal-frontend msigportal-frontend