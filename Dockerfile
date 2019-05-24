FROM debian:9.9-slim

COPY . /opt/MIRACUM-Pipe

RUN apt-get update && apt-get upgrade -y && \
    apt-get install build-essential git-core cmake zlib1g-dev libncurses-dev patch && \
            wget unzip && \
    apt-get remove unzip wget build-essential git-core cmake zlib1g-dev libncurses-dev patch && \
    apt-get autoremove


ENV PATH="/opt/MIRACUM-Pipe/tools/FastQC/:${PATH}"
ENV PATH="/opt/MIRACUM-Pipe/tools/FastQC/:${PATH}"
ENV PATH="/opt/MIRACUM-Pipe/tools/FastQC/:${PATH}"
ENV PATH="/opt/MIRACUM-Pipe/tools/FastQC/:${PATH}"
ENV PATH="/opt/MIRACUM-Pipe/tools/src/FreeC:${PATH}"