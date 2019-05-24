#!/usr/bin/env bash

version_trimmomatic=0.39


###############
# Trimmomatic #
###############

# remove old verion
rm -rf Trimmomatic

# download new version
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${version_trimmomatic}.zip \
    -O trimmomatic.zip

# unzip
unzip -o trimmomatic.zip
rm -f trimmomatic.zip

# rename folder and file (neglect version information)
mv Trimmomatic* Trimmomatic
mv Trimmomatic/trimmomatic-${version_trimmomatic}.jar Trimmomatic/trimmomatic.jar

#########
# FREEC #
#########

# compile FREEC
cd src/FREEC/src && make

# remove object files
rm -f *.o

# copy binary
cd ../
mkdir bin
mv src/freec bin

# go back into tools/
cd ../../


#################
# bam-readcount #
#################

