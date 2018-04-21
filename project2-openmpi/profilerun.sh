#!/bin/bash
set -x

MPIRUN="$(which mpirun)"
$MPIRUN -x VT_NMFILE=pool.nm -x VT_FILE_PREFIX=pool.logs --oversubscribe --hostfile hostfile -np $1 ./pool spec.txt ppmresults/finalbrd
