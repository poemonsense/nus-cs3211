#ifndef MYMPI_H
#define MYMPI_H

#include "mpi.h"

extern MPI_Datatype CompSpecType;
extern MPI_Datatype PoolType;
extern MPI_Datatype LocationType;
extern MPI_Datatype ParticleType;
extern MPI_Datatype RGBPixelType;

int init_mympi();

#define mympi_bcast(buf, count, datatype, root) \
    MPI_Bcast(buf, count, datatype, root, MPI_COMM_WORLD);

#define mympi_isend(buf, count, datatype, dest, tag, request) \
    MPI_Isend(buf, count, datatype, dest, tag, MPI_COMM_WORLD, request)

#define mympi_irecv(buf, count, datatype, source, tag, request) \
    MPI_Irecv(buf, count, datatype, source, tag, MPI_COMM_WORLD, request)

#define mympi_waitall(count, array_of_requests, array_of_statuses) \
    MPI_Waitall(count, array_of_requests, array_of_statuses)

#define mympi_wait(request, status) \
    MPI_Wait(request, status)

#define mympi_send(buf, count, datatype, dest, tag) \
    MPI_Send(buf, count, datatype, dest, tag, MPI_COMM_WORLD)

#define mympi_recv(buf, count, datatype, source, tag, status) \
    MPI_Recv(buf, count, datatype, source, tag, MPI_COMM_WORLD, status)

#define mympi_barrier() \
    MPI_Barrier(MPI_COMM_WORLD);

#endif