#!/usr/bin/env bash

# variables
version_trimmomatic="0.39"
version_picard="2.20.6"
version_VarScan="2.3.9"
version_bedtools="2.28.0"

########
DIR_SCRIPT="$( cd "$(dirname "$0")" ; pwd -P )"

MY_PATH="$DIR_SCRIPT/../"
MY_LD_LIBRARY_PATH=''


# remove old object files
find ${DIR_SCRIPT} -type f -name '*.o' -delete


##########
# FastQC #
##########
cd ${DIR_SCRIPT}/FastQC
ant build
chmod +x bin/fastqc
MY_PATH="$MY_PATH:$DIR_SCRIPT/FastQC/bin/"


####### install from web #######


MY_PATH="$MY_PATH:$DIR_SCRIPT/gatk/"


##########
# picard #
##########
mkdir -p ${DIR_SCRIPT}/picard
cd ${DIR_SCRIPT}/picard

wget https://github.com/broadinstitute/picard/releases/download/${version_picard}/picard.jar \
    -O picard.jar

MY_PATH="$MY_PATH:$DIR_SCRIPT/picard/"


###########
# VarScan #
###########
mkdir -p ${DIR_SCRIPT}/varscan
cd ${DIR_SCRIPT}/varscan
wget https://sourceforge.net/projects/varscan/files/VarScan.v${version_VarScan}.jar \
    -O varscan/VarScan.jar

MY_PATH="$MY_PATH:$DIR_SCRIPT/varscan/"


#############
# bedtools2 #
#############
cd ${DIR_SCRIPT}

wget https://github.com/arq5x/bedtools2/releases/download/v${version_bedtools}/bedtools-${version_bedtools}.tar.gz \
    -O bedtools2.tar.gz

tar -xzf bedtools2.tar.gz
rm -f bedtools2.tar.gz

cd bedtools2 && make

MY_PATH="$MY_PATH:$DIR_SCRIPT/bedtools2/bin/"


##########
# SNPEFF #
##########
cd ${DIR_SCRIPT}
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip -O snpEff.zip

unzip -o snpEff.zip
rm -f snpEff.zip

MY_PATH="$MY_PATH:$DIR_SCRIPT/snpEff/:$DIR_SCRIPT/clinEff/"

###############
# Trimmomatic #
###############
cd ${DIR_SCRIPT}

# download new version
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${version_trimmomatic}.zip \
    -O trimmomatic.zip

# unzip
unzip -o trimmomatic.zip
rm -f trimmomatic.zip

# rename folder and file (neglect version information)
mv Trimmomatic* Trimmomatic
mv Trimmomatic/trimmomatic-${version_trimmomatic}.jar Trimmomatic/trimmomatic.jar

MY_PATH="$MY_PATH:$DIR_SCRIPT/Trimmomatic/"

###############
# Trimmomatic #
###############
cd ${DIR_SCRIPT}



###### COMPILE SUBMODULES #######


#########
# FREEC #
#########
# compile FREEC
cd ${DIR_SCRIPT}/FREEC/src && make

# remove object files
rm -f *.o depend.mk

# copy binary
cd ../
mkdir -p bin
chmod +x src/freec
mv src/freec bin

MY_PATH="$DIR_SCRIPT/FREEC/bin:$MY_PATH"

# add module
cd cd ${DIR_SCRIPT}/FREEC/mappability
wget https://xfer.curie.fr/get/nil/7hZIk1C63h0/hg19_len100bp.tar.gz
tar -xzf hg19_len100bp.tar.gz
rm -f hg19_len100bp.tar.gz


#################
# bam-readcount #
#################
cd ${DIR_SCRIPT}/bam-readcount
mkdir -p build && cd build
cmake ../ && make

# move file
mv ${DIR_SCRIPT}/bam-readcount/build/bin ${DIR_SCRIPT}/bam-readcount/bin
rm -rf ${DIR_SCRIPT}/bam-readcount/build

MY_PATH="$DIR_SCRIPT/bam-readcount/bin:$MY_PATH"

#######
# bwa #
#######
cd ${DIR_SCRIPT}/bwa && make && chmod +x bwa
rm -f *.o

MY_PATH="$DIR_SCRIPT/bwa:$MY_PATH"


##########
# htslib #
##########
cd ${DIR_SCRIPT}/htslib
autoheader     # If using configure, generate the header template...
autoconf       # ...and configure script (or use autoreconf to do both)

./configure    # Optional but recommended, for choosing extra functionality


############
# samtools #
############
cd ${DIR_SCRIPT}/samtools

# create config
autoheader
autoconf -Wno-syntax

# configure and make
./configure

# build samtools and htslib
make

rm -f ${DIR_SCRIPT}/samtools/*.o

MY_PATH="$DIR_SCRIPT/bwa:$MY_PATH"

# htslib is built by samtools
rm -f ${DIR_SCRIPT}/htslib/*.o

# MY_LD_LIBRARY_PATH="$DIR_SCRIPT/htslib:$MY_LD_LIBRARY_PATH"

# add lib folder system wide
echo "$DIR_SCRIPT/htslib" > /etc/ld.so.conf.d/htslib.conf

# write MY_PATH into file
echo "export PATH=$PATH:$MY_PATH" >> /etc/environment
