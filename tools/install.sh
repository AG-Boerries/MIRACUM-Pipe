#!/usr/bin/env bash

# variables
version_trimmomatic=0.39


########
working_dir="$( cd "$(dirname "$0")" ; pwd -P )"

MY_PATH=${working_dir}/../
MY_LD_LIBRARY_PATH=''

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

MY_PATH="$working_dir/Trimmomatic/:$MY_PATH"


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
make all all-htslib

rm -f ${working_dir}/samtools/*.o

MY_PATH="$working_dir/bwa:$MY_PATH"

# htslib is built by samtools
rm -f ${working_dir}/htslib/*.o
MY_LD_LIBRARY_PATH="$working_dir/htslib:$MY_LD_LIBRARY_PATH"




# write MY_PATH into file
echo ${MY_PATH} > ${working_dir}/my_path
echo ${MY_LD_LIBRARY_PATH} > ${working_dir}/my_ld_library_path
