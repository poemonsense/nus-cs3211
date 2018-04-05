#!/bin/bash
set -x

MPIRUN="$(which mpirun)"
$MPIRUN --oversubscribe --hostfile hostfile -np $1 ./pool spec.txt ppmresults/finalbrd
