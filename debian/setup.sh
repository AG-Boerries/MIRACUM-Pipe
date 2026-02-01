#!/usr/bin/env bash

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit
  pwd -P
)

function install_java8()
{
  ## install jre8
  apt-get install -y --no-install-recommends apt-transport-https ca-certificates wget dirmngr gnupg software-properties-common && \
  wget -qO - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add - && \
  add-apt-repository --yes https://packages.adoptium.net/artifactory/deb/ && \
  apt-get update && apt-get install -y --no-install-recommends temurin-8-jdk
}

function install_texlive()
{
  apt-get install -y --no-install-recommends texlive texlive-lang-german texlive-latex-extra
  tlmgr option repository http://mirror.ctan.org/systems/texlive/tlnet #ftp://tug.org/historic/systems/texlive/2018/tlnet-final #http://mirror.ctan.org/systems/texlive/tlnet
  tlmgr install \
    breakurl \
    multirow \
    wrapfig \
    varwidth \
    threeparttable \
    preprint
  tlmgr init-usertree
}

function install_R()
{
  echo "deb http://cloud.r-project.org/bin/linux/debian buster-cran40/" >> /etc/apt/sources.list && \
  apt-key add "/opt/MIRACUM-Pipe/debian/r_key.asc"
  #apt-key adv --keyserver keyserver.ubuntu.com --recv-key B8F25A8A73EACF41
  apt-get update && apt-get install -y --no-install-recommends -t buster-cran40 r-base-dev
  R CMD javareconf
}

# update packages
apt-get update & apt-get upgrade -y

# packages that are required for installation
# libcurl4-openssl-dev
apt-get install -y --no-install-recommends build-essential gcc-multilib libc-dev git-core cmake patch cmake ca-certificates \
  autoconf wget zip unzip zlib1g-dev libbz2-dev liblzma-dev libssl-dev libmariadb-dev tabix libmpfr-dev \
  libncurses5-dev libxml2-dev libcairo2-dev libxt-dev libgit2-dev libcurl4-gnutls-dev libncursesw5-dev libhts-dev libncurses5-dev \
  python python-dev python-numpy python-biopython python-xlrd python-openpyxl libharfbuzz-dev libfribidi-dev libtiff-dev\
  gfortran \
  default-jre \
  ant \
  perl-base \
  python3 python3-pysam python3-pip python3-numpy python3-scipy python3-matplotlib python3-reportlab python3-pandas python3-biopython python3-pyfaidx python3-pyvcf cython python3-setuptools python3-dev libpython3-all-dev python3-future \
  cnvkit \
  libsnappy-java && \
  install_java8 && \
  apt-get -y purge  default-jre default-jdk-headless \
                    openjdk-11-jdk openjdk-11-jdk-headless \
                    openjdk-11-jre openjdk-11-jre-headless && \
  install_R && \
  install_texlive && \
  python3 -m pip install wheel && \
  python3 -m pip install shyaml agfusion==1.4.1 pyensembl==2.2.9 gtfparse==1.3.0 sequenza-utils && \
  apt-get -y autoremove
