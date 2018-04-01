#include <stdio.h>

#include "mympi.h"
#include "spec.h"
#include "pool.h"
#include "mpi.h"


MPI_Datatype CompSpecType;
MPI_Datatype PoolType;

int init_mympi() {
    if (MPI_Type_create_struct(COMPSPEC_COUNT, COMPSPEC_BLOCK_LENGTH,
            COMPSPEC_DISPLACE, COMPSPEC_ELEM_TYPES, &CompSpecType) != MPI_SUCCESS) {
        fprintf(stderr, "[%s:%d] MPI_Type_create_struct failed\n", __FILE__, __LINE__);
        return -1;
    }
    if (MPI_Type_commit(&CompSpecType) != MPI_SUCCESS) {
        fprintf(stderr, "[%s:%d] MPI_Type_commit failed\n", __FILE__, __LINE__);
        return -1;
    }
    if (MPI_Type_create_struct(POOL_COUNT, POOL_BLOCK_LENGTH,
            POOL_DISPLACE, POOL_ELEM_TYPES, &PoolType) != MPI_SUCCESS) {
        fprintf(stderr, "[%s:%d] MPI_Type_create_struct failed\n", __FILE__, __LINE__);
        return -1;
    }
    if (MPI_Type_commit(&PoolType) != MPI_SUCCESS) {
        fprintf(stderr, "[%s:%d] MPI_Type_commit failed\n", __FILE__, __LINE__);
        return -1;
    }
    return 0;
}
