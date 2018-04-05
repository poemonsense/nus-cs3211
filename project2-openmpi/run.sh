#!/bin/bash
set -x

make clean && make "${@:2}"
mpirun --oversubscribe --hostfile hostfile -np $1 ./pool spec.txt ppmresults/finalbrd
make bmp
