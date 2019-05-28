FROM debian:9.9

ADD . /opt/MIRACUM-Pipe

RUN apt-get update && apt-get upgrade -y && \
    /opt/MIRACUM-Pipe/tools/setup.debian.sh && \
    /opt/MIRACUM-Pipe/tools/install.sh && \
    export MY_PATH=`cat /opt/MIRACUM-Pipe/tools/my_path`

ENV PATH="$PATH:$MY_PATH"

RUN echo $PATH