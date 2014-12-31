#!/bin/bash

if [ $# -ne 1 ];then
	echo "Usage: $0 [Release|Debug]"
	exit
fi

if [[ $1 != "Release" && $1 != "Debug" ]];then
	echo "Only Release or Debug allowed, but '$1' is given"
	exit
fi

if [ ! -d $1 ]; then
	mkdir $1
fi

cd $1

if [ -f CMakeCache.txt ]; then
    rm CMakeCache.txt
fi

cmake .. -DCMAKE_INSTALL_PREFIX=~/jlpeng/bin/descriptors -DBUILD_SHARED=OFF -DBABEL_INCLUDE_DIR=~/program/openbabel-2.3.2/include/openbabel-2.0 -DBABEL_LIB_DIR=~/program/openbabel-2.3.2/lib -DCMAKE_BUILD_TYPE=$1
cd ..

