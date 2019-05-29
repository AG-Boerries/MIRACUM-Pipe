FROM debian:9.9

ADD . /opt/MIRACUM-Pipe

RUN apt-get update && apt-get upgrade -y && \
    /opt/MIRACUM-Pipe/tools/setup.debian.sh && \
    /opt/MIRACUM-Pipe/tools/install.sh

ENTRYPOINT ["/opt/MIRACUM-Pipe/interface.sh"]
CMD ["true", "batman", "superman"]

RUN echo $PATH