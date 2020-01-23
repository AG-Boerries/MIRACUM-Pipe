FROM debian:10

ENV HTTP_PROXY "http://proxy.ibsm.uni-freiburg.de:8080"
ENV HTTPS_PROXY "http://proxy.ibsm.uni-freiburg.de:8080"

ADD . /opt/MIRACUM-Pipe

RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils && \
    apt-get upgrade -y && \
    /opt/MIRACUM-Pipe/debian/setup.sh && \
    Rscript /opt/MIRACUM-Pipe/RScripts/install_packages.R && \
    /opt/MIRACUM-Pipe/tools/install.sh

# define volumes
VOLUME /opt/MIRACUM-Pipe/assets/input
VOLUME /opt/MIRACUM-Pipe/assets/output
VOLUME /opt/MIRACUM-Pipe/assets/references

VOLUME /opt/MIRACUM-Pipe/tools/annovar
VOLUME /opt/MIRACUM-Pipe/tools/gatk

VOLUME /opt/MIRACUM-Pipe/databases

VOLUME /opt/MIRACUM-Pipe/conf
