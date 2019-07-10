FROM debian:10

ADD . /opt/MIRACUM-Pipe

RUN apt-get update && apt-get upgrade -y && \
    /opt/MIRACUM-Pipe/tools/setup.debian.sh
    /opt/MIRACUM-Pipe/tools/install.sh && \
    Rscript /opt/MIRACUM-Pipe/RScripts/install_packages.R

ENTRYPOINT ["/opt/MIRACUM-Pipe/interface.sh"]
CMD ["true", "batman", "superman"]

RUN echo $PATH