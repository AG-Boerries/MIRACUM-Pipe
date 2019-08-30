#!/usr/bin/env bash

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit
  pwd -P
)

function install_java8()
{
  ## install jre8
  apt-get install -y apt-transport-https ca-certificates wget dirmngr gnupg software-properties-common && \
  wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | apt-key add - && \
  add-apt-repository --yes https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/ && \
  apt-get update && apt-get install adoptopenjdk-8-openj9
}

function install_texlive()
{
  apt-get install texlive texlive-lang-german
  tlmgr option repository ftp://tug.org/historic/systems/texlive/2018/tlnet-final
  tlmgr install \
    breakurl \
    multirow \
    wrapfig \
    varwidth \
    threeparttable \
    preprint
}
# update packages
apt-get update

# packages that are required for installation
apt-get install -y build-essential gcc-multilib libc-dev git-core cmake patch cmake ca-certificates \
  autoconf wget zip unzip zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev\
  libncurses5-dev libxml2-dev \
  gfortran \
  default-jre \
  ant \
  perl-base \
  r-base-core r-recommended \
  python3 python3-pysam \
  python3-pip \
  libsnappy-java && \
  install_java8 && \
  apt-get -y purge  default-jre default-jdk-headless \
                    openjdk-11-jdk openjdk-11-jdk-headless \
                    openjdk-11-jre openjdk-11-jre-headless && \
  install_texlive && \
  pip3 install shyaml && \
  apt-get purge -y python3-pip && \
  apt-get -y autoremove