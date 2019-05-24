#!/usr/bin/env bash

# packages only required for installation
apt-get remove -y build-essential git-core cmake zlib1g-dev libncurses-dev patch cmake \
 autoconf wget unzip libbz2-dev
apt-get autoremove