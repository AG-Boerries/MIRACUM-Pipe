#!/usr/bin/env bash

DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit
  pwd -P
)

function install_java8()
{
  ## install jdk8
  apt-get install -y apt-transport-https ca-certificates wget dirmngr gnupg software-properties-common && \
  wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | sudo apt-key add - && \
  add-apt-repository --yes https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/ && \
  apt-get update && apt-get install adoptopenjdk-8-openj9 && apt-get remove default-jdk
}


# packages that are required for installation
apt-get install -y build-essential gcc-multilib libc-dev git-core cmake patch cmake  \
  autoconf wget zip unzip zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev\
  libncurses5-dev libxml2-dev \
  gfortran \
  default-jre \
  ant \
  perl-base \
  r-base-core r-recommended \
  texlive \
  python3 python3-pysam \
  python3-pip && \
  pip3 install shyaml && \
  apt-get purge -y python3-pip && \
  install_java8 && \
  apt-get -y purge default-jre && \
  apt-get -y autoremove