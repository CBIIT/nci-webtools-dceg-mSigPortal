# Stage 0 - Build Client
FROM node:latest

COPY client/ /client

WORKDIR /client

RUN npm install \
 && npm run build

# Stage 1 - Copy static assets
FROM centos:latest

RUN dnf -y --setopt=tsflags=nodocs update && \
    dnf -y --setopt=tsflags=nodocs install httpd && \
    dnf clean all

# Simple startup script to avoid some issues observed with container restart
RUN echo -e '#!/bin/bash\nrm -rf /run/httpd/* /tmp/httpd* && exec /usr/sbin/apachectl -DFOREGROUND' > /run-httpd.sh
RUN chmod +x /run-httpd.sh

# Add custom httpd configuration
ADD docker/frontend.conf /etc/httpd/conf.d/frontend.conf

COPY --from=0 /client/build /var/www/html/msigportal

RUN chmod 755 -R /var/www/html/msigportal

EXPOSE 80
EXPOSE 443

CMD ["/run-httpd.sh"]