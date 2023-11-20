#!/usr/bin/bash

if [ ! -f ./lanthanide ]; then
    echo "Binary not found, compiling ..."
    gcc -O3 -o ./lanthanide ../lanthanide.c  -llapack -lblas -lm -lgfortran -lz
fi

if [ ! -f ./vk1k2k3.gz ]; then
    echo "Calculating Vk1k2k3 ..."
    ./lanthanide 0 $1 $2
fi

echo "Calculating energies and eigenvectors ..."
./lanthanide 1 $1 $2