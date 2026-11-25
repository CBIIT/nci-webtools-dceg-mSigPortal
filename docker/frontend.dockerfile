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

FROM public.ecr.aws/amazonlinux/amazonlinux:2023 AS datadog

ARG TARGETARCH

# Install script from https://docs.datadoghq.com/tracing/trace_collection/proxy_setup/httpd/
# mod_datadog is x86_64-only; other arches get an empty dir and httpd skips the module.
RUN mkdir -p /out \
    && if [ "$TARGETARCH" = "amd64" ]; then \
        dnf -y install unzip \
        && cd /tmp \
        && RELEASE_DATA=$(curl -s https://api.github.com/repos/DataDog/httpd-datadog/releases/latest) \
        && DOWNLOAD_URL=$(echo "$RELEASE_DATA" | grep '"browser_download_url".*mod_datadog_artifact.zip"' | sed 's/.*"browser_download_url": *"\([^"]*\)".*/\1/') \
        && curl -Lf -o mod_datadog_artifact.zip "$DOWNLOAD_URL" \
        && unzip -j mod_datadog_artifact.zip -d /out \
        && rm mod_datadog_artifact.zip; \
    fi

FROM public.ecr.aws/amazonlinux/amazonlinux:2023 AS runner

RUN dnf -y update \
    && dnf -y install httpd \
    && dnf clean all

COPY --from=builder /app/client/dist/ /var/www/html/mutational-signatures/

# Add custom httpd configuration
COPY docker/httpd-msigportal.conf /etc/httpd/conf.d/httpd-msigportal.conf
COPY docker/httpd-datadog.conf /etc/httpd/conf.d/00-datadog.conf
COPY --from=datadog /out/ /usr/lib64/httpd/modules/

# forward request and error logs to docker log collector
RUN ln -sf /dev/stdout /var/log/httpd/access_log \
    && ln -sf /dev/stderr /var/log/httpd/error_log

WORKDIR /var/www/html

EXPOSE 80
EXPOSE 443

ENV SERVER_TIMEOUT=900

# Remove Python 3.9 and its dependents (including dnf/yum) from the runtime image.
RUN dnf -y --setopt=protected_packages= --setopt=clean_requirements_on_remove=False remove python3

CMD rm -rf /run/httpd/* /tmp/httpd* \
    && exec /usr/sbin/httpd -DFOREGROUND
