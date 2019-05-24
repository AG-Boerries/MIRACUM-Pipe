FROM debian:9.9-slim

ADD . /opt/MIRACUM-Pipe

RUN apt-get update && apt-get upgrade -y && \
    apt-get install build-essential git-core cmake zlib1g-dev libncurses-dev patch cmake && \
            wget unzip && \
    apt-get remove unzip wget build-essential git-core cmake zlib1g-dev libncurses-dev patch cmake && \
    apt-get autoremove

RUN export MY_PATH=`cat /opt/MIRACUM-Pipe/tools/my_path`
ENV PATH="$MY_PATH:$PATH"

RUN export MY_LD_LIBRARY_PATH=`cat /opt/MIRACUM-Pipe/tools/my_ld_library_path`
ENV LD_LIBRARY_PATH="$MY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"

RUN echo $PATH