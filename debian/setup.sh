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
  apt-get update && apt-get install adoptopenjdk-8-hotspot
}

function install_jdk8()
{
  ## install jre8
   apt-get install -y apt-transport-https ca-certificates wget curl dirmngr gnupg software-properties-common && \
   cp /etc/apt/sources.list /etc/apt/sources.list.bak && \
   echo "deb http://ftp.de.debian.org/debian sid main" >> /etc/apt/sources.list && \
   apt-get update && apt-get -y install openjdk-8-jdk java8-runtime-headless ant && \
   mv /etc/apt/sources.list.bak /etc/apt/sources.list && \
   apt-get update
}

function install_texlive()
{
  apt-get install -y texlive texlive-lang-german texlive-latex-extra
  tlmgr option repository ftp://tug.org/historic/systems/texlive/2018/tlnet-final
  tlmgr install \
    breakurl \
    multirow \
    wrapfig \
    varwidth \
    threeparttable \
    preprint
  tlmgr init-usertree
}
# update packages
apt-get update

# packages that are required for installation
apt-get install -y build-essential gcc-multilib libc-dev git-core cmake patch cmake ca-certificates \
  autoconf wget zip unzip zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libmariadbclient-dev\
  libncurses5-dev libxml2-dev \
  gfortran \
  default-jre \
  ant \
  perl-base \
  r-base-core r-recommended r-cran-latticeextra r-cran-hmisc r-cran-rmysql\
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
  #update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java 
