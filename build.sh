#!/bin/bash

HERE=`pwd`

# rm (possibly conflicting) already installed openbabel
sudo apt remove libopenbabel libopenbabel-dev openbabel

# install open-babel-2.4.0 from sources, system-wide
tar xzvf openbabel-openbabel-2-4-0.tar.gz
cd openbabel-openbabel-2-4-0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr ../
NPROCS=`getconf _NPROCESSORS_ONLN`
make -j $NPROCS
sudo make install

# compile pharao
cd $HERE
mkdir pharao-3.0.4
cd pharao-3.0.4
cmake ../
make -j $NPROCS

# inform about the freshly compiled exe
ls -l pharao-3.0.4/pharao
