#ifndef MYMPI_H
#define MYMPI_H

#include "mpi.h"

extern MPI_Datatype CompSpecType;
extern MPI_Datatype PoolType;

int init_mympi();

#define mympi_bcast(buf, count, datatype, root, comm) \
    MPI_Bcast(buf, count, datatype, root, comm);

#endif