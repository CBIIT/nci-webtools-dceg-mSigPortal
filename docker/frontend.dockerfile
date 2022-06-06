FROM ${FRONTEND_BASE_IMAGE:-quay.io/centos/centos:stream9}

RUN dnf -y update \
    && dnf -y install \
    dnf-plugins-core \
    epel-release \
    && curl -fsSL https://rpm.nodesource.com/setup_16.x | bash - \
    && dnf -y install \
    gcc-c++ \
    httpd \
    nodejs \
    make \
    && dnf clean all

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
    && exec /usr/sbin/httpd -DFOREGROUND

# docker build -t msigportal-frontend -f frontend.dockerfile ~/Projects/msigportal/
# docker run -d -p 8331:80 --name msigportal-frontend msigportal-frontend