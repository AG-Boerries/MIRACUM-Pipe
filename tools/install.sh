#!/usr/bin/env bash

# variables
readonly VERSION_TRIMMOMATIC="0.39"
# 2.20.6
readonly VERSION_PICARD="2.26.4"
readonly VERSION_VARSCAN="2.4.4"
# 2.28.0
readonly VERSION_BEDTOOLS="2.30.0"

########
readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
  pwd -P
)


# remove old object files
find ${DIR_SCRIPT} -type f -name '*.o' -delete


##########
# FastQC #
##########
cd ${DIR_SCRIPT}/FastQC
ant build
chmod +x bin/fastqc


####### install from web #######


##########
# picard #
##########
mkdir -p ${DIR_SCRIPT}/picard
cd ${DIR_SCRIPT}/picard

wget https://github.com/broadinstitute/picard/releases/download/${VERSION_PICARD}/picard.jar \
    -O picard.jar


###########
# VarScan #
###########
mkdir -p ${DIR_SCRIPT}/varscan
cd ${DIR_SCRIPT}/varscan
wget https://github.com/dkoboldt/varscan/raw/master/VarScan.v${VERSION_VARSCAN}.jar \
    -O VarScan.jar


#############
# bedtools2 #
#############
cd ${DIR_SCRIPT}

wget https://github.com/arq5x/bedtools2/releases/download/v${VERSION_BEDTOOLS}/bedtools-${VERSION_BEDTOOLS}.tar.gz \
    -O bedtools2.tar.gz

tar -xzf bedtools2.tar.gz
rm -f bedtools2.tar.gz

cd bedtools2 && make

##########
# SNPEFF #
##########
cd ${DIR_SCRIPT}
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip -O snpEff.zip

unzip -o snpEff.zip
rm -f snpEff.zip

# download database
cd ${DIR_SCRIPT}/snpEff
wget https://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCh37.75.zip -O GRCh37.75.zip
unzip -o GRCh37.75.zip
rm -f GRCh37.37.zip

###############
# Trimmomatic #
###############
cd ${DIR_SCRIPT}

# download new version
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${VERSION_TRIMMOMATIC}.zip \
    -O trimmomatic.zip

# unzip
unzip -o trimmomatic.zip
rm -f trimmomatic.zip

# rename folder and file (neglect version information)
mv Trimmomatic* Trimmomatic
mv Trimmomatic/trimmomatic-${VERSION_TRIMMOMATIC}.jar Trimmomatic/trimmomatic.jar


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

# add module
cd ${DIR_SCRIPT}/FREEC/mappability
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

#######
# bwa-mem2 #
#######
cd ${DIR_SCRIPT}
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 -O bwa-mem2.tar.bz2
tar -xf bwa-mem2.tar.bz2
rm bwa-mem2*.tar.bz2
mv bwa-mem2* bwa-mem2

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

# htslib is built by samtools
rm -f ${DIR_SCRIPT}/htslib/*.o


# add lib folder system wide
echo "$DIR_SCRIPT/htslib" > /etc/ld.so.conf.d/htslib.conf


#################
# fusioncatcher #
#################
cd ${DIR_SCRIPT}
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py
python bootstrap.py --prefix=${DIR_SCRIPT} -t -y

##################
# sequenza-utils #
##################
#pip3 install sequenza-utils

##############
# msisensor2 #
##############
cd ${DIR_SCRIPT}/msisensor2
chmod +x msisensor2

#################
# msisensor-pro #
#################
cd ${DIR_SCRIPT}/msisensor-pro
#wget https://github.com/xjtu-omics/msisensor-pro/raw/master/binary/msisensor-pro
#chmod +x msisensor-pro
./INSTALL

############
# agfusion #
############
cd /opt/MIRACUM-Pipe/databases
agfusion download -g hg38
