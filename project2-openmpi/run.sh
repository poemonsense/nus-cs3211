#!/bin/bash
set -x

MPIRUN="$(which mpirun)"

make clean && make "${@:2}"
$MPIRUN --oversubscribe --hostfile hostfile -np $1 ./pool spec.txt ppmresults/finalbrd
make bmp
