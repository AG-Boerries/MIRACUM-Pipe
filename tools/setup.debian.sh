#!/usr/bin/env bash

# packages that are required for installation
apt-get install -y build-essential gcc-multilib libc-dev git-core cmake patch cmake  \
 autoconf wget zip unzip zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev\
 libncurses5-dev libxml2-dev \
 gfortran \
 default-jre \
 r-base-core r-recommended \
 texlive \
 python3 python3-pysam \
 python3-pip && \
 pip install shyaml && \
 apt-get purge -y python3-pip