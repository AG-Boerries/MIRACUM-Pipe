#!/usr/bin/env bash

# variables
version_trimmomatic=0.39
version_GATK=4.1.2.0


########
working_dir="$( cd "$(dirname "$0")" ; pwd -P )"

MY_PATH="$working_dir/../"
MY_LD_LIBRARY_PATH=''


# remove old object files
find ${working_dir} -type f -name '*.o' -delete


####### install from web #######


########
# GATK #
########
# remove old verion
cd ${working_dir}

# download new version
wget https://github.com/broadinstitute/gatk/releases/download/${version_GATK}/gatk-${version_GATK}.zip \
    -O gatk.zip

# unzip
unzip -o gatk.zip
rm -f gatk.zip

# rename folder and file (neglect version information)
mv gatk* gatk

MY_PATH="$MY_PATH:$working_dir/gatk/"


##########
# picard #
##########
mkdir ${working_dir}/picard
cd ${working_dir}/picard

wget https://github.com/broadinstitute/picard/releases/download/2.20.2/picard.jar \
    -O picard/picard.jar

MY_PATH="$MY_PATH:$working_dir/picard/"


###########
# VarScan #
###########
mkdir -p ${working_dir}/varscan
cd ${working_dir}/varscan
wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar \
    -O varscan/VarScan.jar

MY_PATH="$MY_PATH:$working_dir/varscan/"


#############
# bedtools2 #
#############
cd ${working_dir}

wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz \
    -O bedtools2.tar.gz

tar -xzf bedtools2.tar.gz
rm -f bedtools2.tar.gz

cd bedtools2 && make

MY_PATH="$MY_PATH:$working_dir/bedtools2/bin/"


##########
# SNPEFF #
##########
cd ${working_dir}
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip -O snpEff.zip

unzip -o snpEff.zip
rm -f snpEff.zip

MY_PATH="$MY_PATH:$working_dir/snpEff/:$working_dir/clinEff/"

###############
# Trimmomatic #
###############
cd ${working_dir}

# download new version
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${version_trimmomatic}.zip \
    -O trimmomatic.zip

# unzip
unzip -o trimmomatic.zip
rm -f trimmomatic.zip

# rename folder and file (neglect version information)
mv Trimmomatic* Trimmomatic
mv Trimmomatic/trimmomatic-${version_trimmomatic}.jar Trimmomatic/trimmomatic.jar

MY_PATH="$MY_PATH:$working_dir/Trimmomatic/"

###############
# Trimmomatic #
###############
cd ${working_dir}



###### COMPILE SUBMODULES #######


#########
# FREEC #
#########
# compile FREEC
cd ${working_dir}/FREEC/src && make

# remove object files
rm -f *.o depend.mk

# copy binary
cd ../
mkdir bin
chmod +x src/freec
mv src/freec bin

MY_PATH="$working_dir/FREEC/bin:$MY_PATH"

# add module
cd cd ${working_dir}/FREEC/mappability
wget https://xfer.curie.fr/get/nil/7hZIk1C63h0/hg19_len100bp.tar.gz
tar -xzf hg19_len100bp.tar.gz
rm -f hg19_len100bp.tar.gz


#################
# bam-readcount #
#################
cd ${working_dir}/bam-readcount
mkdir -p build && cd build
cmake ../ && make

# move file
mv ${working_dir}/bam-readcount/build/bin ${working_dir}/bam-readcount/bin
rm -rf ${working_dir}/bam-readcount/build

MY_PATH="$working_dir/bam-readcount/bin:$MY_PATH"

#######
# bwa #
#######
cd ${working_dir}/bwa && make && chmod +x bwa
rm -f *.o

MY_PATH="$working_dir/bwa:$MY_PATH"


##########
# htslib #
##########
cd ${working_dir}/htslib
autoheader     # If using configure, generate the header template...
autoconf       # ...and configure script (or use autoreconf to do both)

./configure    # Optional but recommended, for choosing extra functionality


############
# samtools #
############
cd ${working_dir}/samtools

# create config
autoheader
autoconf -Wno-syntax

# configure and make
./configure

# build samtools and htslib
make

rm -f ${working_dir}/samtools/*.o

MY_PATH="$working_dir/bwa:$MY_PATH"

# htslib is built by samtools
rm -f ${working_dir}/htslib/*.o

# MY_LD_LIBRARY_PATH="$working_dir/htslib:$MY_LD_LIBRARY_PATH"

# add lib folder system wide
echo "$working_dir/htslib" > /etc/ld.so.conf.d/htslib.conf

# write MY_PATH into file
echo "export PATH=$PATH:$MY_PATH" >> /etc/environment
