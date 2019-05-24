#!/usr/bin/env bash

# packages that are required for installation
apt-get install build-essential git-core cmake zlib1g-dev libncurses-dev patch \
        wget unzip


# clean files that were only required for installation
apt-get remove unzip wget build-essential git-core cmake zlib1g-dev libncurses-dev patch
apt-get autoremove
